#![allow(dead_code)]

//extern crate rhai;
//extern crate rand;

pub mod vec;
pub mod color;
pub mod canvas;
pub mod transform;
pub mod shape;
pub mod material;
pub mod camera;
pub mod lua;


use vec::{Point, Ray, Vector};
use canvas::Canvas;
use shape::{Shape, Sphere, Plane, World};
use camera::Camera;
use color::Color;
use transform::Matrix;
use material::{Light, Material, Pattern, 
               StripePattern, GradientPattern, 
               RingPattern, CheckerPattern,
               GridPattern};

use std::f64::consts::PI;
use std::path::Path;

//use rhai::{Engine, RegisterFn};
use rand::Rng;
use simplelog::*;

pub fn floats_equal(f1:f64, f2:f64) -> bool {
    return (f1 - f2).abs() < 0.0001;
}

fn test1() {
    let origin = Point::new(0., 0., -5.);
    let wall_z = 10.;
    let wall_size = 7.;
    let pixel_size = wall_size / 100.;
    let half = wall_size / 2.;
    let mut canvas = Canvas::new(100, 100);
    let color = Color::RED;
    let shape = Sphere::new_with_transform(Matrix::identity().scaling(0.5, 0.5, 0.5));
    
    for i in 0..100 {
        let world_y = half - pixel_size*(i as f64);
        for j in 0..100 {
            let world_x = -half + pixel_size*(j as f64);
            let position = Point::new(world_x, world_y, wall_z);
            let ray = Ray::new(origin, (position.sub(origin)).normalize());
            let intersections = shape.intersect(&ray);
            if intersections.len() > 0 {
                canvas.write_pixel(j, i, color);
            }
        }
    }

    canvas.write_to_file("test1.jpg");
}

fn test2() {
    let origin = Point::new(0., 0., -5.);
    let wall_z = 10.;
    let wall_size = 7.;
    let pixels_per_side = 300;
    let pixel_size = wall_size / pixels_per_side as f64;
    let half = wall_size / 2.;
    let mut canvas = Canvas::new(pixels_per_side, pixels_per_side);
    let shape = Sphere::new_with_transform_and_material
        (Matrix::identity().scaling(1.3, 1.3, 1.3), Material::solid_with_defaults(Color::new(1., 0.2, 1.)));
    let light = Light::new(Color::WHITE, Point::new(-10., 10., -10.));
    for i in 0..pixels_per_side {
        let world_y = half - pixel_size*(i as f64);
        for j in 0..pixels_per_side {
            let world_x = -half + pixel_size*(j as f64);
            let position = Point::new(world_x, world_y, wall_z);
            let ray = Ray::new(origin, (position.sub(origin)).normalize());
            match shape.intersect(&ray).get_hit() {
                Some(hit) => {
                    let hit_point = ray.position(hit.t);
                    let normal = hit.shape.normal_at(&hit_point);
                    let eye = ray.direction.negate();
                    let color = hit.shape.get_material().lighting(&light, Some(&shape), &hit_point, &eye, &normal, false); 
                    canvas.write_pixel(j, i, color);
                }
                None => {} 
            }
        }
    }

    canvas.write_to_file("test2.jpg");
}

fn test3() {
    let mut world = World::new(Light::new(Color::WHITE, Point::new(0., 10., -10.)));
    let shape = Sphere::new_with_transform_and_material(Matrix::identity().scaling(2., 2., 2.).translation(0., 0., 0.),
                                                        Material { color:Some(Color::new(0.1, 1.0, 0.5)), 
                                                                   diffuse:0.7,
                                                                   specular:0.3,
                                                                   ..Material::DEFAULT
                                                        });
    world.add_shape(Box::new(shape));
    
    let camera = Camera::new_with_transform(300, 200, PI/3.0,
                                            Matrix::make_view_transform(Point::new(0., 0., -5.),
                                                                        Point::new(0., 0., 0.),
                                                                        Vector::new(0., 1., 0.)));
    let canvas = camera.render(&world);
    canvas.write_to_file("test3.jpg");
}

/*
* rhai -- not used
fn test4() {

    let mut engine = Engine::new();
    engine.register_type::<Vector>();
    engine.register_type::<Point>();
    engine.register_fn("vector", Vector::new);
    engine.register_fn("point", Point::new);
    engine.register_type::<Matrix>();
    engine.register_fn("identity", Matrix::identity);
    
    println!("{:?}", engine.eval::<Vector>("let v = vector(1.0, 2.0, 3.0); v"));
    println!("{:?}", engine.eval::<Matrix>("let m = identity(); m"));
    
}
*/

/*
fn test5() {

    let mut engine = Engine::new();
    engine.register_type::<Light>();
    engine.register_type::<Sphere>();
    engine.register_type::<Material>();

    engine.register_fn("default_light", <material::Light as Default>::default);
    engine.register_fn("sphere", Sphere::new);
    engine.register_fn("identity", Matrix::identity);
    engine.register_fn("color", Color::color_from_spec);
    engine.register_fn("material", Material::solid);
    
}
*/

fn test6() {
    let mut world = World::new(Default::default());

    let mut pattern = RingPattern::new(Color::GREEN, Color::WHITE);
    pattern.set_transform(Matrix::identity().scaling(0.5, 0.5, 0.5));

    let shape = Sphere::new_with_transform_and_material(Matrix::identity().translation(-0.5, 1., 0.5),
                                                        Material { pattern:Some(Box::new(pattern)),
                                                                  diffuse:0.7,
                                                                  specular:0.3,
                                                                  ..Material::DEFAULT
                                                       });
    world.add_shape(Box::new(shape));

    let mut pattern = GradientPattern::new(Color::RED, Color::BLUE);
    pattern.set_transform(Matrix::identity().scaling(1., 1., 1.));
    let shape = Sphere::new_with_transform_and_material(Matrix::identity().scaling(0.5, 0.5, 0.5).translation(1., 0.7, -4.5),
                                                        Material { pattern:Some(Box::new(pattern)),
                                                                   diffuse:0.7,
                                                                   specular:0.3,
                                                                   ..Material::DEFAULT
                                                        });
    
    world.add_shape(Box::new(shape));

    let mut pattern = StripePattern::new(Color::WHITE, Color::BLACK);
    pattern.set_transform(Matrix::identity().scaling(0.2, 0.2, 0.2));
    
    let shape = Sphere::new_with_transform_and_material(Matrix::identity()
                                                        .scaling(0.8, 0.8, 0.8)
                                                        .translation(-2.5, 0.53, -0.75).
                                                        rotation_x(PI/4.0),
                                                        Material { pattern:Some(Box::new(pattern)),
                                                                   diffuse:0.7,
                                                                   specular:0.3,
                                                                   ..Material::DEFAULT
                                                        });
    
    world.add_shape(Box::new(shape));

    let mut pattern = CheckerPattern::new(Color::WHITE, Color::BLUE);
    pattern.set_transform(Matrix::identity().scaling(1., 1., 1.));
    let plane = Plane::new_with_transform_and_material(Matrix::identity(), 
                                                       Material::pattern_with_defaults(Box::new(pattern))); 
    world.add_shape(Box::new(plane));

    let camera = Camera::new_with_transform(800, 600, PI/2.0,
                                            Matrix::make_view_transform(Point::new(0., 0.5, -5.),
                                                                        Point::new(0., 1., 0.),
                                                                        Vector::new(0., 1., 0.)));
    let canvas = camera.render(&world);
    canvas.write_to_file("test6.jpg");
}

fn test7() {
    let mut world = World::new(Default::default());

    let shape = Sphere::new_with_transform_and_material(Matrix::identity().translation(-0.5, 1., 0.5),
                                                        Material { color:Some(Color::RED),
                                                                   diffuse:0.7,
                                                                   specular:0.3,
                                                                   shininess:1.0,
                                                                  ..Material::DEFAULT
                                                       });
    world.add_shape(Box::new(shape));

    let shape = Sphere::new_with_transform_and_material(Matrix::identity().scaling(0.5, 0.5, 0.5).translation(1., 0.7, -3.5),
                                                        Material { color:Some(Color::GREEN),
                                                                   diffuse:0.7,
                                                                   specular:0.3,
                                                                   shininess:0.5,
                                                                   ..Material::DEFAULT
                                                        });
    
    world.add_shape(Box::new(shape));

    let shape = Sphere::new_with_transform_and_material(Matrix::identity()
                                                        .scaling(0.8, 0.8, 0.8)
                                                        .translation(-2.5, 0.53, -0.75).
                                                        rotation_x(PI/4.0),
                                                        Material { color:Some(Color::BLUE),
                                                                   diffuse:0.7,
                                                                   specular:0.3,
                                                                   shininess:0.2,
                                                                   ..Material::DEFAULT
                                                        });
    
    world.add_shape(Box::new(shape));

    let plane = Plane::new_with_transform_and_material(Matrix::identity(), 
                                                       Material { color:Some(Color::WHITE),
                                                                  diffuse:0.7,
                                                                  specular:0.3,
                                                                  shininess:1.0,
                                                                  ..Material::DEFAULT
                                                        }); 
    world.add_shape(Box::new(plane));

    let camera = Camera::new_with_transform(800, 600, PI/2.0,
                                            Matrix::make_view_transform(Point::new(0., 0.5, -5.),
                                                                        Point::new(0., 1., 0.),
                                                                        Vector::new(0., 1., 0.)));
    let canvas = camera.render(&world);
    canvas.write_to_file("test7.jpg");
}

fn test8() {
    let mut world = World::new(Default::default());

    let shape = Sphere::new_with_transform_and_material(Matrix::identity().translation(-0.5, 0., 0.5),
                                                        Material { color:Some(Color::RED),
                                                                   diffuse:0.1,
                                                                   transparency:1.0,
                                                                   reflectiveness:0.9,
                                                                   refractive_index:1.8,
                                                                   specular:0.9,
                                                                   ambient:0.1,
                                                                  ..Material::DEFAULT
                                                       });
    world.add_shape(Box::new(shape));

    let shape = Sphere::new_with_transform_and_material(Matrix::identity().scaling(0.5, 0.5, 0.5).translation(1., 0.7, -1.5),
                                                        Material { color:Some(Color::GREEN),
                                                                   diffuse:0.1,
                                                                   transparency:1.0,
                                                                   reflectiveness:1.0,
                                                                   ambient:0.1,
                                                                   refractive_index:1.5,
                                                                   specular:0.9,
                                                                   shininess:400.,
                                                                   ..Material::DEFAULT
                                                        });
    
    world.add_shape(Box::new(shape));

    let shape = Sphere::new_with_transform_and_material(Matrix::identity()
                                                        .scaling(0.8, 0.8, 0.8)
                                                        .translation(-1.5, 1., 1.2),
                                                        Material { color:Some(Color::BLUE),
                                                                   diffuse:0.07,
                                                                   transparency:1.0,
                                                                   reflectiveness:0.5,
                                                                   specular:0.3,
                                                                   refractive_index:1.25,
                                                                   ..Material::DEFAULT
                                                        });
    
    world.add_shape(Box::new(shape));

    let mut pattern = GridPattern::new(Color::WHITE, Color::BLACK);
    pattern.set_transform(Matrix::identity().scaling(2.0, 2.0, 2.0));

    let plane = Plane::new_with_transform_and_material(Matrix::identity().translation(0., -1.5, 0.), 
                                                       Material { pattern:Some(Box::new(pattern)),
                                                                  diffuse:0.2,
                                                                  ambient:0.6,
                                                                  reflectiveness:0.5,
                                                                  specular:0.3,
                                                                  ..Material::DEFAULT
                                                        }); 
    world.add_shape(Box::new(plane));

    let mut camera = Camera::new_with_transform(800, 600, 1.3,
                                            Matrix::make_view_transform(Point::new(0., 7., -10.),
                                                                        Point::new(0., 0., 0.),
                                                                        Vector::new(0., 1., 0.)));
    camera.antialiasing_samples = 20;
    let canvas = camera.render_async(&world);
    canvas.write_to_file("test8.jpg");
}

fn test9() {

    let mut world = World::new(Default::default());
    let mut rng = rand::thread_rng();
    
    for _ in 0..30 {
        let mut material = Material::solid_with_defaults(Color::new(rng.gen(), rng.gen(), rng.gen()));
        material.ambient = rng.gen();
        material.diffuse = rng.gen();
        material.specular = 0.0; // rng.gen();
        material.shininess = 0.0; // rng.gen::<f64>()*200.0;
        material.reflectiveness = 0.0; // = rng.gen();
        
        let sphere = Sphere::new_with_transform_and_material
            (Matrix::identity().scaling(0.5, 0.5, 0.5).translation((rng.gen::<f64>()*20.)-10.0, 0., (rng.gen::<f64>()*30.)-15.0),
             material);
        world.add_shape(Box::new(sphere));
    }
    
    world.add_shape(Box::new(Sphere::new_with_transform_and_material
                             (Matrix::identity(),
                              Material::solid_with_defaults(Color::RED))));

    // let mut checked_pattern = CheckerPattern::new(Color::WHITE, Color::BLACK);
    // checked_pattern.set_transform(Matrix::identity().scaling(2.0, 2.0, 2.0));

    let mut pattern = GridPattern::new(Color::color_from_spec("#7983b0").unwrap(), Color::WHITE);
    pattern.set_transform(Matrix::identity().scaling(1.0, 1.0, 1.0));
    
    let plane = Plane::new_with_transform_and_material(Matrix::identity().translation(0., -0.5, 0.), 
                                                       Material { pattern:Some(Box::new(pattern)),
                                                                  diffuse:0.2,
                                                                  ambient:0.6,
                                                                  reflectiveness:0.0,
                                                                  specular:0.0,
                                                                  ..Material::DEFAULT
                                                        }); 
    world.add_shape(Box::new(plane));
    
    let camera = Camera::new_with_transform(800, 600, 0.7,
                                                Matrix::make_view_transform(Point::new(0., 2., -3.),
                                                                            Point::new(0., 0., 0.),
                                                                            Vector::new(0., 1., 0.)));
    //camera.antialiasing_samples = 20;
    let canvas = camera.render_async(&world);
    canvas.write_to_file("test9.jpg");
    
}

fn main() {
    SimpleLogger::init(LevelFilter::Trace, Config::default());
    match lua::render_lua(Path::new("ex2.lua")) {
        Ok(s) => println!("{}", s),
        Err(e) => panic!(e)
    }
}
