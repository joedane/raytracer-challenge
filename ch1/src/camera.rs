


use super::transform::Matrix;
use super::vec::{Point, Ray};
use super::canvas::Canvas;
use super::shape::World;
use std::fmt::Debug;
use rayon::prelude::*;

use itertools::Itertools;

#[derive(Debug)]
pub struct Camera {
    pub hsize:u32,
    pub vsize:u32,
    pub fov:f64,
    pub half_height:f64,
    pub half_width:f64,
    pub pixel_size:f64,
    view_transform:Matrix
}

impl Camera {

    pub const MAX_REFLECTIONS:u8 = 5;

    pub fn new_with_transform(hsize:u32, vsize:u32, fov:f64, m:Matrix) -> Camera {
        let mut c = Camera::new(hsize, vsize, fov);
        c.view_transform = m;
        return c;
    }

    pub fn new(hsize:u32, vsize:u32, fov:f64) -> Camera {
        /*
         * we assume that the screen is one unit in front of the eye, which means that
         * tan(fov/2) is one-half the screen width in world units.
         */
        let half_view = (fov/2.0).tan();
        let aspect = hsize as f64 / vsize as f64;
        let half_width;
        let half_height;
        if aspect < 1.0 {
            half_height = half_view;
            half_width = half_view * aspect;
        } else {
            half_width = half_view;
            half_height = half_view / aspect;
        }
        let pixel_size = (half_width*2.) / hsize as f64;
        Camera { hsize, vsize, fov, half_height, half_width, pixel_size, view_transform:Matrix::identity() }
    }

    pub fn ray_for_pixel(&self, x:u32, y:u32) -> Ray {
        let xoffset = (x as f64 + 0.5) * self.pixel_size;
        let yoffset = (y as f64 + 0.5) * self.pixel_size;
        let world_x = self.half_width - xoffset;
        let world_y = self.half_height - yoffset;

        let vt = &self.view_transform.inverse();
        let pixel = (*vt).transform_point(&Point::new(world_x, world_y, -1.));
        let origin = (*vt).transform_point(&Point::new(0., 0., 0.));
        let direction = pixel.sub(origin).normalize();
        return Ray::new(origin, direction);

    }

    pub fn render(&self, w:&World) -> Canvas {
        let mut canvas = Canvas::new(self.hsize, self.vsize);

        for y in 0..self.vsize-1 {
            for x in 0..self.hsize-1 {
                let ray = self.ray_for_pixel(x, y);
                let color = w.color_at(&ray, Camera::MAX_REFLECTIONS);
                canvas.write_pixel(x, y, color);
            }
        }
        return canvas;
    }

    pub fn render_async(&self, world:&World) -> Canvas {
        let mut canvas = Canvas::new(self.hsize, self.vsize);
        let pixels:Vec<(u32,u32)> = (0..self.hsize).cartesian_product(0..self.vsize).collect();
        let mut results = vec![];
        
        pixels.par_iter().map_with(world, |w, p| {
            let ray = self.ray_for_pixel(p.0, p.1);        
            let color = w.color_at(&ray, Camera::MAX_REFLECTIONS);
            (p, color)
        }).collect_into_vec(&mut results);

        results.iter().for_each(|r| {
            let (pixel, color) = r;
            canvas.write_pixel(pixel.0, pixel.1, *color);
        });
                         
        return canvas;
    }
}


#[cfg(test)]
mod tests {

    use crate::vec::{Vector, Point};
    use crate::shape::Sphere;
    use super::*;
    use crate::*;
    use std::f64::consts::PI;

    #[test]
    fn test_pixel_size1() {
        let c = Camera::new(200, 125, PI/2.0);
        assert!(floats_equal(c.pixel_size, 0.01));
    }

    #[test]
    fn test_pixel_size2() {
        let c = Camera::new(125, 200, PI/2.0);
        println!("pixel size: {:?}", c.pixel_size);
        assert!(floats_equal(c.pixel_size, 0.01));
    }

    #[test]
    fn test_camera1() {
        let c = Camera::new(201,101, PI/2.0);
        let r = c.ray_for_pixel(100, 50);
        assert!(r.origin.approximately_equal(&Point::new(0., 0., 0.)));
        assert!(r.direction.approximately_equal(&Vector::new(0., 0., -1.)));
    }

    #[test]
    fn test_camera2() {
        let c = Camera::new(201,101, PI/2.0);
        let r = c.ray_for_pixel(0, 0);
        assert!(r.origin.approximately_equal(&Point::new(0., 0., 0.)));
        assert!(r.direction.approximately_equal(&Vector::new(0.66519, 0.33259, -0.66851)));
    }


    #[test]
    fn test_camera3() {
        let c = Camera::new_with_transform(201,101, PI/2.0, 
                                           Matrix::identity().translation(0., -2., 5.).rotation_y(PI/4.0));
        
        let r = c.ray_for_pixel(100, 50);
        assert!(r.origin.approximately_equal(&Point::new(0., 2., -5.)));
        let n = 2.0_f64.sqrt()/2.0;
        assert!(r.direction.approximately_equal(&Vector::new(n, 0., -n)));
    }

    #[test]
    fn test_render1() {
        let world:World = Default::default();
        let from = Point::new(0., 0., -5.);
        let to = Point::new(0., 0., 0.);
        let up = Vector::new(0., 1., 0.);
        let camera = Camera::new_with_transform(11, 11, PI/2., 
                                                Matrix::make_view_transform(from, to, up));
        let c = camera.render(&world);
        let test_color = c.get_pixel(5, 5);
        println!("color: {:?}", test_color);
        assert!(floats_equal(test_color.red, 0.38066));
        assert!(floats_equal(test_color.green, 0.47583));
        assert!(floats_equal(test_color.blue, 0.2855));
        
    }

    #[test]
    fn test_render_jd1() {
        let mut world = World::new(Default::default());
        let shape = Sphere::new_with_transform(1.0, Matrix::identity().scaling(2.0, 2.0, 2.0).translation(0., 1.01, 0.));
        world.add_shape(Box::new(shape));
        
        let from = Point::new(0., 0., -5.);
        let to = Point::new(0., 0., 0.);
        let up = Vector::new(0., 1., 0.);
        let camera = Camera::new_with_transform(101, 101, PI/2., 
                                                Matrix::make_view_transform(from, to, up));
        let c = camera.render(&world);
        let test_color = c.get_pixel(50, 50);
        println!("color: {:?}", test_color);
        c.write_to_file("test.jpg");
    }

    #[test]
    fn test_async1() {
        let c = Camera::new(20, 10, 1.5);
        let w:World = Default::default();
        c.render_async(&w);
    }
}
