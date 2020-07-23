
use array2d::Array2D;
use num::Num;
use num::traits::Zero;
use super::vec::{Vector, Point};


pub fn multiply_arrays<T: Num + Copy>(a1:&Array2D<T>, a2:&Array2D<T>) -> Array2D<T> {
    let mut new_m = Array2D::filled_with(Zero::zero(), 4, 4);
    for row in 0..4 {
        for col in 0..4 {
            let mut sum = Zero::zero();
            for i in 0..4 {
                let v = a1[(row, i)] * a2[(i, col)];
                sum = sum + v;
            }
            new_m[(row, col)] = sum;
        }
    }
    new_m
}

#[derive(Debug, Clone)]
pub struct Matrix {
    matrix:Array2D<f64>,
    determinant:Option<f64>
}

impl Matrix {
        
    pub fn multiply(&self, other:&Matrix) -> Matrix {
        Matrix::new(multiply_arrays(&self.matrix, &other.matrix))  
    }

    pub fn new(m:Array2D<f64>) -> Self {
        let d = Matrix::determinant(&m);
        Matrix { matrix:m, determinant: if d.abs() > Vector::EPSILON { Some(d) } else { None }}
    }

    pub fn new_from_rows(rows: &[std::vec::Vec<f64>]) -> Self {
        Matrix::new(Array2D::from_rows(rows))
    }

    pub fn identity() -> Matrix {
        let m = Array2D::from_row_major(&vec!
            [1., 0., 0., 0., 
             0., 1., 0., 0.,
             0., 0., 1., 0.,
             0., 0., 0., 1.], 4, 4);
        Matrix::new(m)
    }

    pub fn translation(&self, x:f64, y:f64, z:f64) -> Matrix {
        let m = Array2D::from_row_major(&vec!
            [1., 0., 0., x, 
             0., 1., 0., y,
             0., 0., 1., z,
             0., 0., 0., 1.], 4, 4);
        Matrix::new(m).multiply(self)
    }

    pub fn scaling(&self, x:f64, y:f64, z:f64) -> Matrix {
        let m = Array2D::from_row_major(&vec!
            [x, 0., 0., 0., 
             0., y, 0., 0.,
             0., 0., z, 0.,
             0., 0., 0., 1.], 4, 4);
        Matrix::new(m).multiply(self)
    }

    pub fn rotation_x(&self, r:f64) -> Matrix {
        let m = Array2D::from_row_major(&vec!
            [1., 0., 0., 0., 
             0., r.cos(), -r.sin(), 0.,
             0., r.sin(), r.cos(), 0.,
             0., 0., 0., 1.], 4, 4);
        Matrix::new(m).multiply(self)
    }

    pub fn rotation_y(&self, r:f64) -> Matrix {
        let m = Array2D::from_row_major(&vec!
            [r.cos(), 0., r.sin(), 0., 
             0., 1., 0., 0.,
            -r.sin(), 0., r.cos(), 0.,
             0., 0., 0., 1.], 4, 4);
        Matrix::new(m).multiply(self)
    }

    pub fn rotation_z(&self, r:f64) -> Matrix {
        let m = Array2D::from_row_major(&vec!
            [r.cos(), -r.sin(), 0., 0., 
             r.sin(), r.cos(), 0., 0.,
             0., 0., 1., 0.,
             0., 0., 0., 1.], 4, 4);
        Matrix::new(m).multiply(self)
    }

    pub fn shearing(&self, xy:f64, xz:f64, yx:f64, yz:f64, zx:f64, zy:f64) -> Matrix {
        let m = Array2D::from_row_major(&vec!
            [1., xy, xz, 0., 
             yx, 1., yz, 0.,
             zx, zy, 1., 0.,
             0., 0., 0., 1.], 4, 4);
        Matrix::new(m).multiply(self)
    }

    pub fn transform_vector(&self, p:&Vector) -> Vector {
        /*
        return Vector::new(self.matrix.row_iter(0).zip(v).fold(0.0, |acc, (l, r)| acc + l*r),
        self.matrix.row_iter(1).zip(v).fold(0.0, |acc, (l, r)| acc + l*r),
        self.matrix.row_iter(2).zip(v).fold(0.0, |acc, (l, r)| acc + l*r));

    */
    return Vector::new(
        self.matrix[(0,0)]*p.x + self.matrix[(0, 1)]*p.y + self.matrix[(0, 2)]*p.z,
        self.matrix[(1,0)]*p.x + self.matrix[(1, 1)]*p.y + self.matrix[(1, 2)]*p.z,
        self.matrix[(2, 0)]*p.x + self.matrix[(2, 1)]*p.y + self.matrix[(2, 2)]*p.z        
    );

    }

    pub fn transform_point(&self, p:&Point) -> Point {
        return Point::new(
            self.matrix[(0,0)]*p.x + self.matrix[(0, 1)]*p.y + self.matrix[(0, 2)]*p.z + self.matrix[(0, 3)],
            self.matrix[(1,0)]*p.x + self.matrix[(1, 1)]*p.y + self.matrix[(1, 2)]*p.z + self.matrix[(1, 3)],
            self.matrix[(2, 0)]*p.x + self.matrix[(2, 1)]*p.y + self.matrix[(2, 2)]*p.z + self.matrix[(2, 3)]
            );
     }

     pub fn determinant(m:&Array2D<f64>) -> f64 {
         if m.num_columns() != m.num_rows() {
             panic!("Can only take a determinant of a square matrix");
         }
        if m.num_columns() == 2 {
            return m[(0,0)] * m[(1,1)] - m[(0, 1)] * m[(1, 0)]
        } else {
            let mut det:f64 = 0.0;
            for col in 0..m.num_columns() {
                det += m[(0, col)] * Matrix::cofactor(m, 0, col);
            }
            return det;
        }
     }

     pub fn submatrix(m:&Array2D<f64>, i:usize, j:usize) -> Array2D<f64> {
        if m.num_columns() <= 2 {
            panic!("Matrix is too small"); 
        }
        let new_size = m.num_columns() - 1;
        let mut new_m:Array2D<f64> = Array2D::filled_with(0.0, new_size, new_size);
        for row in 0..new_size {
            for col in 0..new_size {
                new_m[(row, col)] = m[(if row < i {row} else {row+1}, if col < j {col} else {col+1})];
            }
        }
        new_m
     }

     pub fn minor(m:&Array2D<f64>, i:usize, j:usize) -> f64 {
         Matrix::determinant(&Matrix::submatrix(m, i, j))
     }

     pub fn cofactor(m:&Array2D<f64>, i:usize, j:usize) -> f64 {
        let mut minor = Matrix::minor(m, i, j);
        if (i+j) % 2 == 1 {
            minor = -1.0*minor;
        }
        return minor;
     }

     pub fn is_invertable(&self) -> bool {
        self.determinant.is_some()
     }

    pub fn inverse(&self) -> Matrix {
        match self.determinant {
            None => panic!("Matrix is not invertable"),
            Some(d) => {
                let l = self.matrix.num_columns();
                let mut new_m = Array2D::filled_with(0.0, l, l);
                for row in 0..l {
                    for col in 0..l {
                        let c = Matrix::cofactor(&self.matrix, row, col);
                        new_m[(col, row)] = c / d;
                    }
                }
                Matrix::new(new_m)
            }
        }
    }

    pub fn transpose(&self) -> Matrix {
        let mut m:Array2D<f64> = Array2D::filled_with(0.0, 
                                            self.matrix.num_columns(),
                                            self.matrix.num_rows());
        for row in 0..m.num_rows() {
            for col in 0..m.num_columns() {
                m[(row, col)] = self.matrix[(col, row)];
            }
        }
        Matrix::new(m) 
    }

    pub fn make_view_transform(from:Point, to:Point, up:Vector) -> Matrix {
        let forward = to.sub(from).normalize();
        let upn = up.normalize();
        let left = forward.cross(&upn);
        let true_up = left.cross(&forward);
        
        let m = Array2D::from_row_major(&vec!
            [left.x, left.y, left.z, 0.,
             true_up.x, true_up.y, true_up.z, 0.,
            -(forward.x), -(forward.y), -(forward.z), 0.,
            0., 0., 0., 1.], 4, 4);

            return Matrix::new(multiply_arrays(&m, &Matrix::identity().translation(-from.x, -from.y, -from.z).matrix));
    }
}

impl std::cmp::PartialEq for Matrix {

    fn eq(&self, other: &Self) -> bool {
        self.matrix == other.matrix
    }

}

#[cfg(test)]
mod tests {

    use crate::floats_equal;
    use crate::vec::{Vector, Point};
    use super::{multiply_arrays, Matrix}; 
    use array2d::Array2D;


    #[test]
    fn test_make_array() {
        let m1 = Array2D::from_rows(
            &vec![vec![1, 2, 3, 4],
                            vec![2, 5, 89, 1],
                            vec![43, 1, -10, 100]]);

        assert_eq!(m1[(0, 0)], 1);
        assert_eq!(m1[(2, 0)], 43);
    }

    #[test]
    fn test_array_eq() {
        let m1 = Array2D::from_rows(
            &vec![vec![1, 2, 3],
                            vec![4, 5, 6],
                            vec![7, 8, 9]]);

        let m2 = Array2D::from_rows(
            &vec![vec![1, 2, 3],
                            vec![4, 5, 6],
                            vec![7, 8, 9]]);

        let m3 = Array2D::from_rows(
            &vec![vec![1, 2, 3],
                            vec![4, 5, 6],
                            vec![7, 8, 10]]);

        assert!(m1 == m2);
        assert!(m1 != m3);
    }

    #[test]
    fn test_array_mul() {

        let m1 = Array2D::from_rows(
            &vec![vec![1, 2, 3, 4],
                            vec![5, 6, 7, 8],
                            vec![9, 8, 7, 6],
                            vec![5, 4, 3, 2]]);

        let m2 = Array2D::from_rows(
            &vec![vec![-2, 1, 2, 3],
                            vec![3, 2, 1, -1],
                            vec![4, 3, 6, 5],
                            vec![1, 2, 7, 8]]);

        let m3 = Array2D::from_rows(
            &vec![vec![20, 22, 50, 48],
                            vec![44, 54, 114, 108],
                            vec![40, 58, 110, 102],
                            vec![16, 26, 46, 42]]);


        assert!(multiply_arrays(&m1, &m2) == m3);

    }
/*
    #[test]
    fn test_array_vec_mul() {
        let m = Array2D::from_rows(
            &vec![vec![1.0, 2.0, 3.0],
                            vec![2.0, 4.0, 4.0],
                            vec![8.0, 6.0, 4.0]]);


        let v = Vector::new(1.0, 2.0, 3.0);
        
        assert!(v.matrix_mul(Matrix::new(m)) == Vector::new(14.0, 22.0, 32.0));
        
    }
    */

    #[test]
    fn test_array_transpose() {

        let m = Matrix::new_from_rows(&vec![
                vec![1.0, 2.0, 3.0],
                vec![4.0, 5.0, 6.0],
                vec![7.0, 8.0, 9.0]]);

        let mt = Matrix::new_from_rows(&vec![
                vec![1.0, 4.0, 7.0],
                vec![2.0, 5.0, 8.0],
                vec![3.0, 6.0, 9.0]]);

        assert!(m.transpose() == mt);
    }

    #[test]

    fn test_array_determinant1() {
        let m = Matrix::new_from_rows(&vec![
                        vec![1.0, 2.0, 6.0],
                        vec![-5.0, 8.0, -4.0],
                        vec![2.0, 6.0, 4.0]]);

        let det:f64 = m.determinant.unwrap();
        assert!(floats_equal(det, -196.0));
    }

        #[test]
        fn test_array_determinant2() {
        
        }

        #[test]
    fn test_array_invert() {
        let m = Matrix::new_from_rows(&vec![
                        vec![-5.0, 2.0, 6.0, -8.0],
                        vec![1.0, -5.0, 1.0, 8.0],
                        vec![7.0, 7.0, -6.0, -7.0],
                        vec![1.0, -3.0, 7.0, 4.0]]);

        let inverse : Array2D<&str> = Array2D::from_rows(&vec![
                vec!["0.21805", "0.45113", "0.24060", "-0.04511"],
                vec!["-0.80827", "-1.45677", "-0.44361", "0.52068"],
                vec!["-0.07895", "-0.22368", "-0.05263", "0.19737"],
                vec!["-0.52256", "-0.81391", "-0.30075", "0.30639"]]);

        assert!(m.inverse().
                matrix.elements_row_major_iter().zip
                    (inverse.elements_row_major_iter()).
                        all(|(val, check)| format!("{:.5}", val) == *check));
        
    }


    #[test]
    fn test_transform() {
        let m = Matrix::identity()
            .rotation_x(std::f64::consts::PI/2.)
            .scaling(5., 5., 5.)
            .translation(10., 5., 7.);
        let p = Point::new(1., 0., 1.);
        let tp = m.transform_point(&p);
        assert!(floats_equal(tp.x, 15.));
        assert!(floats_equal(tp.y, 0.));
        assert!(floats_equal(tp.z, 7.));


    }

    #[test]
    fn test_view_transform1() {
        let m = Matrix::make_view_transform(Point::new(0., 0., 0.), 
                                            Point::new(0., 0., -1.),
                                            Vector::new(0., 1., 0.));
        assert!(m.matrix == Matrix::identity().matrix);
    }

    #[test]
    fn test_view_transform2() {
        let m = Matrix::make_view_transform(Point::new(0., 0., 0.), 
                                            Point::new(0., 0., 1.),
                                            Vector::new(0., 1., 0.));
        assert!(m.matrix == Matrix::identity().scaling(-1., 1., -1.).matrix);
    }

    #[test]
    fn test_view_transform3() {
        let m = Matrix::make_view_transform(Point::new(0., 0., 8.), 
                                            Point::new(0., 0., 0.),
                                            Vector::new(0., 1., 0.));
        println!("view transform: {:?}", m.matrix);
        println!("compare: {:?}", Matrix::identity().translation(0., 0., -8.).matrix);
        assert!(m.matrix == Matrix::identity().translation(0., 0., -8.).matrix);
    }

    #[test]
    fn test_view_transform4() {
        let _m = Matrix::make_view_transform(Point::new(1., 3., 2.), 
                                             Point::new(4., -2., 8.),
                                             Vector::new(1., 1., 0.));
        let _compare: Array2D<&str> = Array2D::from_rows(&vec![
            vec!["-0.50709", "0.50709", "0.67612", "-2.36643"],
            vec!["0.76772", "0.60609", "0.12122", "-2.82843"],
            vec!["-0.35857", "0.59761", "-0.71714", "0.00000"],
            vec!["0.00000", "0.00000", "0.00000", "1.00000"]]);
      
            assert!(_m.matrix.elements_row_major_iter().zip
                (_compare.elements_row_major_iter()).
                    all(|(val, check)| format!("{:.5}", val) == *check));
      
            /*
         * this test appears to be passing by examination of the output.  Need a better
        * way to compare matrices.
        
        println!("{:?}", m.matrix);
        println!("{:}", m.matrix.mapv(|f| format!("{:.5}", f)));
        assert!(m.matrix.mapv(|f| format!("{:.5}", f)) == compare);
         */
    }
   
}
