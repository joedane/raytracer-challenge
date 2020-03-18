#![allow(dead_code)]

pub mod vec;
pub mod color;
pub mod canvas;
pub mod transform;
pub mod shape;
pub mod material;
pub mod camera;

use vec::{Point, Ray, Vector};
use canvas::Canvas;
use shape::{Shape, Sphere, World};
use camera::Camera;
use color::Color;
use transform::Matrix;
use material::{Light, Material};

use std::f32::consts::PI;

pub fn floats_equal(f1:f32, f2:f32) -> bool {
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
    let shape = Sphere::new_with_transform(1., Matrix::identity().scaling(1., 0.5, 1.));
    
    for i in 0..100 {
        let world_y = half - pixel_size*(i as f32);
        for j in 0..100 {
            let world_x = -half + pixel_size*(j as f32);
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
    let pixel_size = wall_size / pixels_per_side as f32;
    let half = wall_size / 2.;
    let mut canvas = Canvas::new(pixels_per_side, pixels_per_side);
    let shape = Sphere::new_with_transform_and_material
        (1., Matrix::identity(), Material::default_with_color(Color::new(1., 0.2, 1.)));
    let light = Light::new(Color::WHITE, Point::new(-10., 10., -10.));
    for i in 0..pixels_per_side {
        let world_y = half - pixel_size*(i as f32);
        for j in 0..pixels_per_side {
            let world_x = -half + pixel_size*(j as f32);
            let position = Point::new(world_x, world_y, wall_z);
            let ray = Ray::new(origin, (position.sub(origin)).normalize());
            match shape.intersect(&ray).get_hit() {
                Some(hit) => {
                    let hit_point = ray.position(hit.t);
                    let normal = hit.shape.normal_at(&hit_point);
                    let eye = ray.direction.negate();
                    let color = hit.shape.get_material().lighting(&light, &hit_point, &eye, &normal); 
                    canvas.write_pixel(j, i, color);
                }
                None => {} 
            }
        }
    }

    canvas.write_to_file("test2.jpg");
}

fn test3() {
    let mut world = World::new(Default::default());
    world.add_shape(Box::new
                    (Sphere {
                        r:1.0, 
                        transform:Matrix::identity().translation(-0.5, 1., 0.5),
                        material:Material { color:Color::new(0.1, 1.0, 0.5), 
                                            diffuse:0.7,
                                            specular:0.3,
                                            ..Material::DEFAULT
                                            }
                    }));
    let camera = Camera::new_with_transform(300, 200, PI/3.0,
                                            Matrix::make_view_transform(Point::new(0., 1.5, -10.),
                                                                        Point::new(0., 0., 0.),
                                                                        Vector::new(0., 1., 0.)));
    let canvas = camera.render(&world);
    canvas.write_to_file("test3.jpg");
}

fn test4() {

}

fn main() {
    test3();
}
