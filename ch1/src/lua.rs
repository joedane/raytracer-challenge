
use std::fs;
use std::fs::File;
use std::path::Path;
use image::gif::Encoder;
use rlua::prelude::*;
use rlua::Error::FromLuaConversionError;
use rlua::{UserData, UserDataMethods};

use crate::shape::{World, Sphere, Plane};
use crate::material::{Light, Material, Pattern, CheckerPattern, StripePattern, GridPattern};
use crate::transform::Matrix;
use crate::camera::Camera;
use crate::color::Color;
use crate::vec::{Point, Vector};

pub struct GifEncoder {
    encoder: Encoder<std::fs::File>,
}

impl GifEncoder {

    fn new(f: std::fs::File) -> GifEncoder {
        GifEncoder {
            encoder: Encoder::new(f)
        }
    }
}
 
impl UserData for GifEncoder {

    fn add_methods<'lua, M: UserDataMethods<'lua, Self>>(methods: &mut M) {

        methods.add_method_mut("AddFrame", |_, this, (worldtable, cameratable): (LuaTable, LuaTable)| {

            let world = world_from_table(worldtable)?;
            let camera = camera_from_table(cameratable)?;
            let canvas = camera.render_async(&world);
            canvas.frame_to_file(&mut this.encoder);
            Ok(())
        });

        methods.add_method("Finish", |_, _this, ()| {
            Ok(())
        });

    }
}

pub fn render_lua(script: &Path) -> LuaResult<&str> {

    let lua = Lua::new();
    
    let _result = lua.context::<_, LuaResult<&str>>(|lua_ctx| {
        let globals = lua_ctx.globals();

        let renderfn = lua_ctx.create_function(|_, (worldtable, cameratable, 
                                                    outfile): (LuaTable, LuaTable, String)| {

            println!("FOO!");
            let world = world_from_table(worldtable)?;
            println!("bar!");
            let camera = camera_from_table(cameratable)?;

            let canvas = camera.render_async(&world);

            canvas.write_to_file(&outfile);

            Ok("Inside")
        }).unwrap();
        globals.set("Render", renderfn)?;

        lua_ctx.scope(|scope| {

            globals.set("StartAnimation",
                        scope.create_function(|_, outfile:String| {
                            Ok(GifEncoder::new(File::create(outfile).unwrap()))
                        }).unwrap()
            ).unwrap();

            let scriptdata = &fs::read(script).unwrap();
            lua_ctx.load(scriptdata).set_name(script.to_str().unwrap()).unwrap().exec().unwrap();

        });
        Ok("Not sure what to use here")
    });
    Ok("Also not sure about this")
}


fn point_from_table(pointtable: &LuaTable) -> LuaResult<Point> {
    Ok(Point::new(pointtable.get("x")?,
                  pointtable.get("y")?,
                  pointtable.get("z")?))
}

fn vector_from_table(pointtable: &LuaTable) -> LuaResult<Vector> {
    Ok(Vector::new(pointtable.get("x")?,
                   pointtable.get("y")?,
                   pointtable.get("z")?))
}

fn color_from_table(colortable: &LuaTable) -> LuaResult<Color> {
    Ok(Color::new(colortable.get("r")?,
                  colortable.get("g")?,
                  colortable.get("b")?))
}

fn pattern_from_table(patterntable: &LuaTable) -> LuaResult<Box<dyn Pattern>> {

    let m:Matrix = transform_from_table(patterntable)?; 

    match patterntable.get::<_, String>("type")?.as_str() {
        "checks" => {
            let color_a = color_from_table(&patterntable.get("color_a")?)?;
            let color_b = color_from_table(&patterntable.get("color_b")?)?;
            let mut p = CheckerPattern::new(color_a, color_b);
            p.set_transform(m);
            return Ok(Box::new(p))
            },
        "stripes" => {
            let color_a = color_from_table(&patterntable.get("color_a")?)?;
            let color_b = color_from_table(&patterntable.get("color_b")?)?;
            let mut p = StripePattern::new(color_a, color_b);
            p.set_transform(m);
            return Ok(Box::new(p))
        },
        "grid" => {
            let mut p = GridPattern::new(Color::WHITE, Color::BLACK);
            p.set_transform(m);
            return Ok(Box::new(p));
        },
        patterntype => {
            return Err(FromLuaConversionError {
                from: "pattern table",
                to: "pattern",
                message: Some(format!("invalid pattern type: {}", patterntype))
            })
        }
    }
}

fn light_from_table(lighttable: &LuaTable) -> LuaResult<Light> {
    Ok(Light::new(color_from_table(&lighttable.get("color")?)?,
                  point_from_table(&lighttable.get("position")?)?))
}

fn lights_from_table(lightstable: &LuaTable) -> LuaResult<Light> {
    light_from_table(&lightstable.get(1)?)
}

fn float_value(v: LuaValue) -> LuaResult<f64> {
    match v {
        LuaValue::Number(n) => Ok(n),
        LuaValue::Integer(i) => Ok(i as f64),
        _ => Err(LuaError::RuntimeError(format!("Invalid number: {:?}", v)))
    }
}

fn u32_value(v: LuaValue) -> LuaResult<u32> {
    match v {
        LuaValue::Integer(n) => {
            if n >= 0 && n <= std::i32::MAX as i64 {
                Ok(n as u32)
            } else {
                Err(LuaError::RuntimeError(format!("Number out of bounds: {}", n)))
            }
        },
        _ => Err(LuaError::RuntimeError(format!("Invalid number: {:?}", v)))
    }
}

fn u8_value(v: LuaValue) -> LuaResult<u8> {
    match v {
        LuaValue::Integer(n) => {
            if n >= 0 && n <= std::u8::MAX as i64 {
                Ok(n as u8)
            } else {
                Err(LuaError::RuntimeError(format!("Number out of bounds: {}", n)))
            }
        },
        _ => Err(LuaError::RuntimeError(format!("Invalid number: {:?}", v)))
    }
}

fn material_from_table(table: &LuaTable) -> LuaResult<Material> {
    let mut material:Material = Default::default();
    if table.contains_key("material")? {

        match table.get::<_, LuaTable>("material") {
            Ok(materialtable) => {
                for pair in materialtable.pairs::<String, LuaValue>() {
                    let (key, value) = pair?;
                    match key.as_str() {
                        "ambient" => { material.ambient = float_value(value)? }, 
                        "diffuse" => { material.diffuse = float_value(value)? },
                        "specular" => { material.specular = float_value(value)? },
                        "shininess" => {  material.shininess = float_value(value)? },
                        "reflectiveness" => {  material.reflectiveness = float_value(value)?},
                        "transparency" => { material.transparency = float_value(value)? },
                        "refractive_index" => { material.refractive_index = float_value(value)? },
                        "color" => { 
                            if let LuaValue::Table(t) = value { 
                                material.color = Some(color_from_table(&t)?) 
                            } else {
                                return Err(LuaError::RuntimeError(format!("invalid color")));
                            }
                        },
                        "pattern" => {
                            if let LuaValue::Table(t) = value {
                                material.pattern = Some(pattern_from_table(&t)?);
                            } else {
                                return Err(LuaError::RuntimeError(format!("invali pattern")));
                            }
                        },
                        _ => { return Err(LuaError::RuntimeError(format!("Invalid material property: {}", key))); } 
                    }
                }
            }
            Err(e) => { return Err(e); } 
        }
    }
    if let Ok(colortable) =  table.get::<_, LuaTable>("color") {
        if let Ok(color) = color_from_table(&colortable) {
            material.color = Some(color);
        }
    } else if let Ok(patterntable) = table.get::<_, LuaTable>("pattern") {
        if let Ok(pattern) = pattern_from_table(&patterntable) {
            material.pattern = Some(pattern);
        }
    }
    Ok(material)
}

fn camera_from_table(table: LuaTable) -> LuaResult<Camera> {
    
    let position = point_from_table(&table.get("position")?)?;
    let lookat = point_from_table(&table.get("lookat")?)?;
    let up = vector_from_table(&table.get("up")?)?;
    let mut camera = Camera::new_with_transform(u32_value(table.get("screenwidth")?)?,
                                                u32_value(table.get("screenheight")?)?,
                                                float_value(table.get("fov")?)?,
                                                Matrix::make_view_transform(position, lookat, up)); 
    if table.contains_key("samples")? {
        camera.set_samples(u8_value(table.get("samples")?)?);
    }
    Ok(camera)
}

fn transform_from_table(table: &LuaTable) -> LuaResult<Matrix> {
    let mut transform = Matrix::identity();
   
    if table.contains_key("rotate_x")? {
        let r = table.get("rotate_x")?;
        let s = float_value(r)?;
        transform = transform.rotation_x(s);
    }

    if table.contains_key("rotate_y")? {
        let r = table.get("rotate_y")?;
        let s = float_value(r)?;
        transform = transform.rotation_y(s);
    } 

    if table.contains_key("rotate_z")? {
        let r = table.get("rotate_z")?;
        let s = float_value(r)?;
        transform = transform.rotation_z(s);
    }

    if table.contains_key("scale")? {
        let scale = table.get("scale")?;
        let s = float_value(scale)?;
        transform = transform.scaling(s, s, s);
    }

    let positiontable = table.get("position")?;
    if let Some(t) = positiontable {
        let origin = point_from_table(&t)?;
        transform = transform.translation(origin.x, origin.y, origin.z);
    }
    Ok(transform)
}

fn sphere_from_table(shapetable: &LuaTable) -> LuaResult<Sphere> {
    let material = material_from_table(shapetable)?;
    let transform = transform_from_table(shapetable)?; 
    Ok(Sphere::new_with_transform_and_material(transform, material))
}

fn plane_from_table(shapetable: &LuaTable) -> LuaResult<Plane> {
    let material = material_from_table(shapetable)?;
    let transform = transform_from_table(shapetable)?;
    Ok(Plane::new_with_transform_and_material(transform, material))
}


fn world_from_table(worldtable: LuaTable) -> LuaResult<World> {
    let lights = lights_from_table(&worldtable.get("lights")?)?;
    let mut world = World::new(lights);
    let shapestable: LuaTable = worldtable.get("shapes")?;
    for shapetable in shapestable.sequence_values::<LuaTable>() {
        let unwrapped_shapetable = shapetable.unwrap(); // FIXME
        match unwrapped_shapetable.get::<&str, String>("type")?.as_str() {
            "sphere" => {
                world.add_shape(Box::new(sphere_from_table(&unwrapped_shapetable)?));
                println!("added sphere");
            },
            "plane" => {
                world.add_shape(Box::new(plane_from_table(&unwrapped_shapetable)?));
                println!("added plane");
            },
            shapetype => return Err(FromLuaConversionError {
                from: "shape table",
                to: "Shape",
                message: Some(format!("Invalid shape type: {}", shapetype))
            }),
        }
    } 
    Ok(world)
}
