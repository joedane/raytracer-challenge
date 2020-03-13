
use std::ops::{Add, Sub, Mul, Neg};

extern crate ndarray;

use ndarray::prelude::*;
use crate::transform::Matrix;


#[derive(Debug, PartialEq)]
pub struct Vector {
    pub x: f32,
    pub y: f32,
    pub z : f32
}

impl Vector {

    pub fn new(x:f32, y:f32, z:f32) -> Vector {
        Vector {
            x:x, y:y, z:z
        }

    }

    pub fn from_array(a:Array1<f32>) -> Vector {
        Vector::new(a[0], a[1], a[2])
    }

    pub fn add(&self, other:Vector) -> Vector {
        Vector {
            x:self.x + other.x,
            y:self.y + other.y,
            z:self.z + other.z
        }
    }

    pub fn sub(&self, other:Vector) -> Vector {
        Vector {
            x:self.x - other.x,
            y:self.y - other.y,
            z:self.z - other.z
        }
    }

    pub fn negate(&self) -> Vector {
        Vector {
            x:-self.x,
            y:-self.y,
            z:-self.z
        }
    }

    pub fn mul(&self, m:f32) -> Vector {
        Vector {
            x:self.x * m,
            y:self.y * m,
            z:self.z * m
        }
    }

    pub fn div(&self, d:f32) -> Vector {
        Vector {
            x:self.x / d,
            y:self.y / d,
            z:self.z / d
        }
    }

    pub fn magnitude(&self) -> f32 {
        (self.x*self.x + self.y*self.y + self.z*self.z).sqrt()
    }

    pub fn normalize(&self) -> Vector {
        let m = self.magnitude();
        Vector {
            x:self.x / m,
            y:self.y / m,
            z:self.z / m
        }
    }

    pub fn dot(&self, other:&Vector) -> f32 {
        self.x * other.x +
            self.y * other.y +
            self.z * other.z
    }

    pub fn cross(&self, other:&Vector) -> Vector {
        Vector {
            x:self.y * other.z - self.z * other.y,
            y:self.z * other.x - self.x * other.z,
            z:self.x * other.y - self.y * other.x
        }
    }

    pub fn matrix_mul(&self, matrix:&Array2<f32>) -> Vector {
        let v = Array::from(vec!(self.x, self.y, self.z));
        let result = matrix.dot(&v);
        Vector {
            x:result[0],
            y:result[1],
            z:result[2]
        }
    }

    pub fn transform(&self, m:&Matrix) -> Vector {
        return m.transform_vector(self);
    }

}

#[derive(Debug, Copy, PartialEq)]
pub struct Point {
    pub x: f32,
    pub y: f32,
    pub z: f32
}

impl Point {

    pub fn new(x:f32, y:f32, z:f32) -> Point {
        Point {
            x:x, y:y, z:z
        }
    }

    pub fn add(&self, other:Point) -> Vector {
        Vector::new(self.x + other.x,
                    self.y + other.y,
                    self.z + other.z)
    }

    fn add_vec(&self, other:Vector) -> Point {
        Point::new(self.x + other.x,
                   self.y + other.y,
                   self.z + other.z)
    }

    pub fn sub(&self, other:Point) -> Vector {
        Vector::new(self.x - other.x,
                    self.y - other.y,
                    self.z - other.z)
    }

    fn sub_vec(&self, other:Vector) -> Point {
        Point::new(self.x - other.x,
                   self.y - other.y,
                   self.z - other.z)
    }

    pub fn from_array(a:Array1<f32>) -> Point {
        Point::new(a[0], a[1], a[2])
    }

    pub fn transform(&self, m:&Matrix) -> Point {
        return m.transform_point(self);
    }

}

#[derive(Debug)]
pub struct Ray {
    pub origin: Point,
    pub direction: Vector
}

impl Ray {

    pub fn new(o: Point, d: Vector) -> Ray {
        Ray {
            origin:o,
            direction:d
        }
    }

    pub fn position(&self, t:f32) -> Point {
        return self.origin.add_vec(self.direction.mul(t));
    }

    pub fn transform(&self, m:&Matrix) -> Ray {
        return Ray::new(m.transform_point(&self.origin),
                        m.transform_vector(&self.direction));
    }
    
}

#[derive(Debug)]
struct GenericVector<T> {
    x: T, 
    y: T,
    z: T
}

#[derive(Debug)]
struct GenericPoint<T> {
    x: T, 
    y: T,
    z: T
}


impl<T: Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Neg<Output = T> + Copy> GenericVector<T> {
    
    fn add(&self, other:GenericVector<T>) -> GenericVector<T> {
        GenericVector {
            x:self.x + other.x, 
            y:self.y + other.y,
            z:self.z + other.z
        }
    } 
    
    fn sub(&self, other:GenericVector<T>) -> GenericVector<T> {
        GenericVector {
            x:self.x - other.x, 
            y:self.y - other.y,
            z:self.z - other.z
        }
    } 

    fn negate(&self) -> GenericVector<T> {
        GenericVector {
            x:-self.x,
            y:-self.y,
            z:-self.z
        }
    }
}


impl<T: Add<Output = T> + Sub<Output = T> + Copy> GenericPoint<T> {
    
    fn sub(&self, other:GenericPoint<T>) -> GenericVector<T> {
        GenericVector {x:self.x - other.x,
                       y:self.y - other.y,
                       z:self.z - other.z
        }
    }

    fn sub_vec(&self, other:GenericVector<T>) -> GenericPoint<T> {
        GenericPoint {
            x:self.x - other.x,
            y:self.y - other.y,
            z:self.z - other.z
        }   
    }

}


#[cfg(test)]
mod tests {

    //use super::*;
    use crate::floats_equal;
    use crate::transform::Matrix;
    use super::Vector;
    use super::Ray;
    use super::Point;

    #[test]
    fn test_make_point() {
        let p = Point{x:1.0, y:2.0, z:3.0};        
        assert_eq!(p.x, 1.0)
    }

    #[test]
    fn test_negate_vec() {
        let v = Vector {x:1.0, y:2.0, z:3.0};
        let nv = v.negate();
        assert_eq!(nv.x, -1.0);
    }

    #[test]
    fn test_add_vec() {
        let v1 = Vector{x:1.0, y:2.0, z:3.0};
        let v2 = Vector{x:10.0, y:20.0, z:30.0};
        let v3 = v1.add(v2);
        assert_eq!(v3.x, 11.0);
    }

    #[test]
    fn test_add_f32_vec() {
        let v1 = Vector{x:1.0, y:2.0, z:3.0};
        let v2 = Vector{x:10.0, y:20.0, z:30.0};
        let v3 = v1.add(v2);
        assert!(floats_equal(v3.x, 11.0));
    }

    #[test]
    fn test_sub_vec() {
        let v1 = Vector{x:1.0, y:2.0, z:3.0};
        let v2 = Vector{x:10.0, y:20.0, z:30.0};
        let v3 = v1.sub(v2);
        assert_eq!(v3.x, -9.0);
    }

    #[test]
    fn test_magnitude() {
        let v1 = Vector{x:2.0, y:5.0, z:4.0};
        let m = v1.magnitude();
        assert_eq!(m, 6.708204);
    }

    #[test]
    fn test_dot() {
        let v1 = Vector{x:1.0, y:2.0, z:3.0};
        let v2 = Vector{x:2.0, y:3.0, z:4.0};
        assert!(floats_equal(v1.dot(&v2), 20.0));
        assert!(floats_equal(v2.dot(&v1), 20.0));
    }
    
    #[test]
    fn test_cross() {
        let v1 = Vector{x:1.0, y:2.0, z:3.0};
        let v2 = Vector{x:2.0, y:3.0, z:4.0};
        let v3 = v1.cross(&v2);
        let v4 = v2.cross(&v1);
        assert!(floats_equal(v3.x, -1.0));
        assert!(floats_equal(v4.x, 1.0));
    }

    #[test]
    fn test_transform_point1() {
        let p = Point::new(-4., 6., 8.);
        let m = Matrix::identity()
            .scaling(2., 3., 4.);
        let tp = p.transform(&m);

        assert_eq!(tp, Point::new(-8., 18., 32.));
        
    }

    #[test]
    fn test_translate_ray1() {
        let r = Ray::new(Point::new(1., 2., 3.),
                         Vector::new(0., 1., 0.));
        let m = Matrix::identity()
            .translation(3., 4., 5.);
        
        let rt = r.transform(&m);

        assert_eq!(rt.origin, Point::new(4., 6., 8.));
        assert_eq!(rt.direction, Vector::new(0., 1., 0.));

    }
    
    #[test]
    fn test_translate_ray2() {
        let r = Ray::new(Point::new(1., 2., 3.),
                         Vector::new(0., 1., 0.));
        let m = Matrix::identity()
            .scaling(2., 3., 4.);
        
        let rt = r.transform(&m);

        assert_eq!(rt.origin, Point::new(2., 6., 12.));
        assert_eq!(rt.direction, Vector::new(0., 3., 0.));

    }

}