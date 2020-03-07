
extern crate ndarray;

//use ndarray_linalg::solve::{Determinant, Inverse};

use ndarray::prelude::*;
use ndarray_linalg::Inverse;
use super::vec::{Vector, Point};

#[derive(Debug)]
pub struct Matrix {
    matrix:Array2<f32>
}

impl Matrix {
    
    pub fn identity() -> Matrix {
        let m = Array::from_shape_vec((4, 4), vec![1., 0., 0., 0., 
                                                   0., 1., 0., 0.,
                                                   0., 0., 1., 0.,
                                                   0., 0., 0., 1.]).unwrap();
        Matrix { matrix:m }
    }

    pub fn translation(&self, x:f32, y:f32, z:f32) -> Matrix {
        let m = Array::from_shape_vec((4, 4), vec![1., 0., 0., x, 
                                                   0., 1., 0., y,
                                                   0., 0., 1., z,
                                                   0., 0., 0., 1.]).unwrap();
        Matrix { matrix:m.dot(&self.matrix) }
    }

    pub fn scaling(&self, x:f32, y:f32, z:f32) -> Matrix {
        let m = Array::from_shape_vec((4, 4), vec![x, 0., 0., 0., 
                                                   0., y, 0., 0.,
                                                   0., 0., z, 0.,
                                                   0., 0., 0., 1.]).unwrap();
        Matrix { matrix:m.dot(&self.matrix) }
    }

    pub fn rotation_x(&self, r:f32) -> Matrix {
        let m = Array::from_shape_vec((4, 4), vec![1., 0., 0., 0., 
                                                   0., r.cos(), -r.sin(), 0.,
                                                   0., r.sin(), r.cos(), 0.,
                                                   0., 0., 0., 1.]).unwrap();
        Matrix { matrix:m.dot(&self.matrix) }
    }

    pub fn rotation_y(&self, r:f32) -> Matrix {
        let m = Array::from_shape_vec((4, 4), vec![r.cos(), 0., r.sin(), 0., 
                                                   0., 1., 0., 0.,
                                                   -r.sin(), 0., r.cos(), 0.,
                                                   0., 0., 0., 1.]).unwrap();
        Matrix { matrix:m.dot(&self.matrix) }
    }

    pub fn rotation_z(&self, r:f32) -> Matrix {
        let m = Array::from_shape_vec((4, 4), vec![r.cos(), -r.sin(), 0., 0., 
                                                   r.sin(), r.cos(), 0., 0.,
                                                   0., 0., 1., 0.,
                                                   0., 0., 0., 1.]).unwrap();
        Matrix { matrix:m.dot(&self.matrix) }
    }

    pub fn shearing(&self, xy:f32, xz:f32, yx:f32, yz:f32, zx:f32, zy:f32) -> Matrix {
        let m = Array::from_shape_vec((4, 4), vec![1., xy, xz, 0., 
                                                   yx, 1., yz, 0.,
                                                   zx, zy, 1., 0.,
                                                   0., 0., 0., 1.]).unwrap();
        Matrix { matrix:m.dot(&self.matrix) }
    }

    pub fn transform_vector(&self, v:&Vector) -> Vector {
        return Vector::from_array(self.matrix.dot(&arr1(&[v.x, v.y, v.z, 0.])));
    }

    pub fn transform_point(&self, p:&Point) -> Point {
        return Point::from_array(self.matrix.dot(&arr1(&[p.x, p.y, p.z, 1.])));
    }

    pub fn inverse(&self) -> Matrix {
        println!("M: {:?}", self.matrix);
        return Matrix {
            matrix: self.matrix.inv().unwrap()
        };
    }
}


#[cfg(test)]
mod tests {

    use crate::floats_equal;
    use crate::vec::{Vector, Point};
    use super::Matrix;
    use ndarray::prelude::*;
    use ndarray_linalg::{Determinant, Inverse};
 

    #[test]
    fn test_make_array() {
        let m1 = arr2(&[[1, 2, 3, 4],
                        [2, 5, 89, 1],
                        [43, 1, -10, 100]]);

        assert_eq!(m1[[0, 0]], 1);
        assert_eq!(m1[[2, 0]], 43);
    }

    #[test]
    fn test_array_eq() {
        let m1 = arr2(&[[1, 2, 3],
                        [4, 5, 6],
                        [7, 8, 9]]);

        let m2 = arr2(&[[1, 2, 3],
                        [4, 5, 6],
                        [7, 8, 9]]);

        let m3 = arr2(&[[1, 2, 3],
                        [4, 5, 6],
                        [7, 8, 10]]);

        assert!(m1 == m2);
        assert!(m1 != m3);
    }

    #[test]
    fn test_array_mul() {

        let m1 = arr2(&[[1, 2, 3, 4],
                        [5, 6, 7, 8],
                        [9, 8, 7, 6],
                        [5, 4, 3, 2]]);

        let m2 = arr2(&[[-2, 1, 2, 3],
                        [3, 2, 1, -1],
                        [4, 3, 6, 5],
                        [1, 2, 7, 8]]);

        let m3 = arr2(&[[20, 22, 50, 48],
                        [44, 54, 114, 108],
                        [40, 58, 110, 102],
                        [16, 26, 46, 42]]);


        assert!(m1.dot(&m2) == m3);

    }

    #[test]
    fn test_array_vec_mul() {
        let m = arr2(&[[1.0, 2.0, 3.0],
                       [2.0, 4.0, 4.0],
                       [8.0, 6.0, 4.0]]);


        let v = Vector::new(1.0, 2.0, 3.0);
        
        assert!(v.matrix_mul(&m) == Vector::new(14.0, 22.0, 32.0));
        
    }
    
    #[test]
    fn test_array_transpose() {

        let m = arr2(&[[1.0, 2.0, 3.0],
                       [2.0, 3.0, 4.0]]);

        let mt = arr2(&[[1.0, 2.0],
                        [2.0, 3.0],
                        [3.0, 4.0]]);

        assert!(m.t() == mt);
    }

    #[test]
    fn test_array_determinant() {
        let m = arr2(&[[1.0, 2.0, 6.0],
                       [-5.0, 8.0, -4.0],
                       [2.0, 6.0, 4.0]]);

        let det:f32 = m.det().unwrap();
        assert!(floats_equal(det, -196.0));
    }

    #[test]
    fn test_array_invert() {
        let m = arr2(&[[-5.0, 2.0, 6.0, -8.0],
                       [1.0, -5.0, 1.0, 8.0],
                       [7.0, 7.0, -6.0, -7.0],
                       [1.0, -3.0, 7.0, 4.0]]);
        let inverse = arr2(&[["0.21805", "0.45113", "0.24060", "-0.04511"],
                             ["-0.80827", "-1.45677", "-0.44361", "0.52068"],
                             ["-0.07895", "-0.22368", "-0.05263", "0.19737"],
                             ["-0.52256", "-0.81391", "-0.30075", "0.30639"]]);

        assert!(m.inv().unwrap().mapv(|f| format!("{:.5}", f)) == inverse);

    }

    #[test]
    fn test_transform() {
        let m = Matrix::identity()
            .rotation_x(std::f32::consts::PI/2.)
            .scaling(5., 5., 5.)
            .translation(10., 5., 7.);
        let p = Point::new(1., 0., 1.);
        let tp = m.transform_point(&p);
        assert!(floats_equal(tp.x, 15.));
        assert!(floats_equal(tp.y, 0.));
        assert!(floats_equal(tp.z, 7.));


    }


}
