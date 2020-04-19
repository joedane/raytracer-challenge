

use std::ops::{Add, Sub, Mul, Neg};
extern crate ndarray;

use ndarray::prelude::*;
use crate::transform::Matrix;
use crate::floats_equal;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z : f64
}

impl Vector {

    pub const EPSILON: f64 = 0.00000001;
    
    pub fn new(x:f64, y:f64, z:f64) -> Vector {
        Vector {
            x:x, y:y, z:z
        }

    }

    pub fn from_array(a:Array1<f64>) -> Vector {
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

    pub fn mul(&self, m:f64) -> Vector {
        Vector {
            x:self.x * m,
            y:self.y * m,
            z:self.z * m
        }
    }

    pub fn div(&self, d:f64) -> Vector {
        Vector {
            x:self.x / d,
            y:self.y / d,
            z:self.z / d
        }
    }

    pub fn magnitude(&self) -> f64 {
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

    pub fn dot(&self, other:&Vector) -> f64 {
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

    pub fn matrix_mul(&self, matrix:&Array2<f64>) -> Vector {
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

    pub fn reflect(&self, n:&Vector) -> Vector {
        return self.sub(n.mul(2. * self.dot(n)));
    }

    pub fn approximately_equal(&self, v:&Vector) -> bool {
        return floats_equal(self.x, v.x) &&
            floats_equal(self.y, v.y) &&
            floats_equal(self.z, v.z);
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl Point {

    pub const ZERO:Point = Point::new(0., 0., 0.);

    pub const fn new(x:f64, y:f64, z:f64) -> Point {
        Point {
            x:x, y:y, z:z
        }
    }

    pub fn add(&self, other:Point) -> Vector {
        Vector::new(self.x + other.x,
                    self.y + other.y,
                    self.z + other.z)
    }

    pub fn add_vec(&self, other:Vector) -> Point {
        Point::new(self.x + other.x,
                   self.y + other.y,
                   self.z + other.z)
    }

    pub fn sub(&self, other:Point) -> Vector {
        Vector::new(self.x - other.x,
                    self.y - other.y,
                    self.z - other.z)
    }

    pub fn sub_vec(&self, other:Vector) -> Point {
        Point::new(self.x - other.x,
                   self.y - other.y,
                   self.z - other.z)
    }

    pub fn from_array(a:Array1<f64>) -> Point {
        Point::new(a[0], a[1], a[2])
    }

    pub fn transform(&self, m:&Matrix) -> Point {
        return m.transform_point(self);
    }

    pub fn approximately_equal(&self, p:&Point) -> bool {
        return floats_equal(self.x, p.x) &&
            floats_equal(self.y, p.y) &&
            floats_equal(self.z, p.z);
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

    pub fn position(&self, t:f64) -> Point {
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
    fn test_add_f64_vec() {
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
        assert!(floats_equal(m, 6.708204));
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

    #[test]
    fn test_vector_reflect1() {
        let v = Vector::new(1., -1., 0.);
        let n = Vector::new(0., 1., 0.);
        let r = v.reflect(&n);
        assert!(floats_equal(r.x, 1.));
        assert!(floats_equal(r.y, 1.));
        assert!(floats_equal(r.z, 0.));
        
    }

    #[test]
    fn test_vector_reflect2() {
        let v = Vector::new(0., -1., 0.);
        let f = (2f64).sqrt()/2.;
        let n = Vector::new(f, f, 0.);
        let r = v.reflect(&n);
        assert!(floats_equal(r.x, 1.));
        assert!(floats_equal(r.y, 0.));
        assert!(floats_equal(r.z, 0.));
        
    }

}
