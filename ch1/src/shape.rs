

use super::vec::{Ray, Point};
use super::transform::Matrix;

use std::fmt::Debug;
use std::option::Option;
use std::ops::Index;

#[derive(Debug)]
pub struct Intersection<'a> {
    t:f32,
    shape:&'a dyn Shape
}

impl<'a> Intersection<'a> {

    fn new(t:f32, shape:&'a dyn Shape) -> Intersection {
        Intersection {
            t:t,
            shape:shape
        }
    }

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        println!("BAR");
        write!(f, "XXX")
    } 
    
}

impl<'a> PartialEq for Intersection<'a> {

    fn eq(&self, other: &Self) -> bool {
        self.t == other.t
    }
}


pub struct Intersections<'a> {
    _i: Vec<Intersection<'a>>,
    _nearest_positive: Option<usize>
}

impl<'a> Intersections<'a> {

    const EMPTY:Intersections<'a> = Intersections { _i: Vec::new(), _nearest_positive: None };
    
    fn new(i:Intersection<'a>) -> Intersections<'a> {
        let mut v = Vec::new();
        let i_t = i.t;
        v.push(i);
        Intersections {
            _i:v,
            _nearest_positive: if i_t >= 0. { Some(0) } else { None }
        }
    }

    fn add(&mut self, i:Intersection<'a>) {
        if i.t >= 0.  {
            match self._nearest_positive {
                Some(l) => if i.t < self._i[l].t {
                    self._nearest_positive = Some(self._i.len());
                }
                None => {
                    self._nearest_positive = Some(self._i.len());
                }
            }

        }
        self._i.push(i);
    }

    pub fn len(&self) -> usize {
        return self._i.len();
    }

    pub fn get_hit(&self) -> Option<&Intersection> {
        match self._nearest_positive {
            Some(l) => { Some(&self._i[l]) }
            None => None
        }
    }
}

impl<'a> Index<usize> for Intersections<'a> {

    type Output = Intersection<'a>;

    fn index(&self, i: usize) -> &Self::Output {
        return &self._i[i];
    }
}


pub trait Shape : Debug {

    fn intersect(&self, ray: &Ray) -> Intersections;

    //fn get_transform(&self) -> Matrix;

}

#[derive(Debug)]
pub struct Circle {
    r: f32,
    transform: Matrix
}

impl Circle {

    pub fn new(r:f32) -> Circle {
        Circle::new_with_transform(r, Matrix::identity())
    }

    pub fn new_with_transform(r:f32, m:Matrix) -> Circle {
        Circle {
            r:r,
            transform:m.inverse()
        }
    }
}


impl Shape for Circle {

    fn intersect(&self, ray: &Ray) -> Intersections {
        let ray = ray.transform(&self.transform);
        let sphere_to_ray = ray.origin.sub(Point::new(0., 0., 0.));
        let a = ray.direction.dot(&ray.direction);
        let b = 2. * ray.direction.dot(&sphere_to_ray);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.;
        let discriminant = (b*b) - 4. * a * c;
        if discriminant < 0. {
            return Intersections::EMPTY;
        } else { 
            let mut i = Intersections::new(Intersection::new((-b - discriminant.sqrt()) / (2. * a), self));
            i.add(Intersection::new((-b + discriminant.sqrt()) / (2. * a), self));
            return i;
        }
    }
}


#[cfg(test)]
mod tests {

    use crate::vec::{Vector, Ray, Point};
    use super::*;
    use crate::*;

    #[test]
    fn test_intersect1() {
        let c = Circle::new(1.);
        let ray = Ray::new(Point::new(0., 0., -5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);
        assert_eq!(v.len(), 2);
        
        assert!(floats_equal(v[0].t, 4.) || floats_equal(v[0].t, 6.));
        assert!(floats_equal(v[1].t, 4.) || floats_equal(v[1].t, 6.));

    }

    #[test]
    fn test_intersect2() {
        let c = Circle::new(1.);
        let ray = Ray::new(Point::new(0., 1., -5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert!(floats_equal(v[0].t, 5.) && floats_equal(v[0].t, 5.));
        
    }

    #[test]
    fn test_intersect3() {
        let c = Circle::new(1.);
        let ray = Ray::new(Point::new(0., 2., -5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert_eq!(v.len(), 0);
    }

    #[test]
    fn test_intersect4() {
        let c = Circle::new(1.);
        let ray = Ray::new(Point::new(0., 0., 0.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert_eq!(v.len(), 2);
        
        assert!(floats_equal(v[0].t, -1.) || floats_equal(v[0].t, 1.));
        assert!(floats_equal(v[1].t, -1.) || floats_equal(v[1].t, 1.));
    }

    #[test]
    fn test_intersect5() {
        let c = Circle::new(1.);
        let ray = Ray::new(Point::new(0., 0., 5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert_eq!(v.len(), 2);
        
        assert!(floats_equal(v[0].t, -4.) || floats_equal(v[0].t, -6.));
        assert!(floats_equal(v[1].t, -4.) || floats_equal(v[1].t, -6.));
    }

    #[test]
    fn test_hits1() {
        let c = Circle::new(1.);
        let i1 = Intersection::new(1., &c);
        let i2 = Intersection::new(2., &c);
        let mut is = Intersections::new(i1);
        is.add(i2);

        assert_eq!(is.get_hit().unwrap().t, 1.);
    }

    #[test]
    fn test_hits2() {
        let c = Circle::new(1.);
        let i1 = Intersection::new(-1., &c);
        let i2 = Intersection::new(-2., &c);
        let mut is = Intersections::new(i1);
        is.add(i2);

        assert!(is.get_hit().is_none());
    }

    #[test]
    fn test_shape_transform1() {
        let c = Circle::new_with_transform(1., Matrix::identity().scaling(2., 2., 2.));
        let r = Ray::new(Point::new(0., 0., -5.),
                         Vector::new(0., 0., 1.));
        let is = c.intersect(&r);
        assert_eq!(2, is.len());
        
        assert!(floats_equal(is[0].t, 3.) || floats_equal(is[0].t, 7.));
        assert!(floats_equal(is[1].t, 3.) || floats_equal(is[1].t, 7.));
    }

    #[test]
    fn test_shape_transform2() {
        let c = Circle::new_with_transform(1., Matrix::identity().translation(5., 0., 0.));
        let r = Ray::new(Point::new(0., 0., -5.),
                         Vector::new(0., 0., 1.));
        let is = c.intersect(&r);
        assert_eq!(0, is.len());
        
    }

}
