
pub mod vec;
pub mod color;
pub mod canvas;
pub mod transform;
pub mod shape;

use vec::{Point, Ray};
use canvas::Canvas;
use shape::{Shape, Circle};
use color::Color;

fn test1() {
    let origin = Point::new(0., 0., -5.);
    let wall_z = 10.;
    let wall_size = 7.;
    let pixel_size = wall_size / 100.;
    let half = wall_size / 2.;
    let canvas = Canvas::new(100, 100);
    let color = Color::RED;
    let shape = Circle::new(1.);
    
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
}

fn main() {
    test1();
}
