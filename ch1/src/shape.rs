

use super::vec::{Ray, Point, Vector};
use super::transform::Matrix;
use super::material::{Material, Light};
use super::color::Color;


use std::fmt::Debug;
use std::vec::Vec;
use std::option::Option;
use std::ops::Index;
use std::borrow::{Borrow, BorrowMut};
use itertools::Itertools;


pub trait Shape : Debug {

    fn intersect(&self, ray: &Ray) -> Intersections;

    fn normal_at(&self, p: &Point) -> Vector;
    //fn get_transform(&self) -> Matrix;

    fn get_material(&self) -> &Material;

    fn get_material_mut(&mut self) -> &mut Material;
}

pub struct CachedVectors {
    t:f32, 
    point:Point,
    eyev:Vector,
    normal:Vector,
    inside:bool
}

impl CachedVectors {
    
    pub fn new(t:f32, point:Point, eyev:Vector, normal:Vector) -> Self{
        let inside = normal.dot(&eyev) < 0.0;
        
        return CachedVectors { 
            t:t, 
            point:point, 
            eyev:eyev, 
            normal: if inside { normal.negate() } else { normal },
            inside:inside
        };
    }
}

#[derive(Debug)]
pub struct Intersection<'a> {
    pub t:f32,
    pub shape:&'a dyn Shape
}

impl<'a> Intersection<'a> {

    fn new(t:f32, shape:&'a dyn Shape) -> Intersection<'a> {
        Intersection {
            t:t,
            shape:shape
        }
    }

    fn compute_vectors(&self, r:&Ray) -> CachedVectors {
        let p = r.position(self.t);
        let eyev = r.direction.negate();
        let normal = self.shape.normal_at(&p);
        return CachedVectors::new(self.t, p, eyev, normal);
    }

}

impl<'a> PartialEq for Intersection<'a> {

    fn eq(&self, other: &Self) -> bool {
        self.t == other.t
    }
}


pub struct Intersections<'a> {
    _i: Vec<Intersection<'a>>
}

impl<'a> Intersections<'a> {

    const EMPTY:Intersections<'a> = Intersections { _i: Vec::new() };
    
    pub fn new(i:Intersection<'a>) -> Intersections<'a> {
        Intersections {
            _i:vec![i]
        }
    }

    pub fn add(&mut self, i:Intersection<'a>) {
        self._insert_sorted(i);
    }

    fn _insert_sorted(&mut self, x:Intersection<'a>) {
        if self._i.len() == 0 {
            self._i.push(x);
        }
        else {
            for i in 0..self._i.len() {
                if x.t < self._i[i].t {
                    self._i.insert(i, x);
                    return;
                }
            }
            self._i.push(x);
        }
    }

    pub fn len(&self) -> usize {
        return self._i.len();
    }

    pub fn merge(&mut self, is:Intersections<'a>) {
        for i in is._i {
            self.add(i);
        }
    }

    pub fn get_hit(&self) -> Option<&Intersection> {
        if self.len() == 0 {
            return None;
        }
        else {
            for i in 0..self.len() {
                if self._i[i].t >= 0.0 {
                    return Some(&self._i[i]);
                }
            }
            return None;
        }
    }

    pub fn get_intersections_sorted(&self) -> std::vec::IntoIter<&Intersection<'a>> {
        return self._i.iter().sorted_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    }
}

impl<'a> IntoIterator for Intersections<'a> {

    type Item = Intersection<'a>;
    type IntoIter = std::vec::IntoIter<Intersection<'a>>;

    fn into_iter(self) -> std::vec::IntoIter<Intersection<'a>> {
        return self._i.into_iter();
    }
}

impl<'s, 'a> IntoIterator for &'s Intersections<'a> {

    type Item = &'s Intersection<'a>;
    type IntoIter = std::slice::Iter<'s, Intersection<'a>>;

    fn into_iter(self) -> std::slice::Iter<'s, Intersection<'a>> {
        return self._i.iter();
    }
    
}

impl<'a> Index<usize> for Intersections<'a> {

    type Output = Intersection<'a>;

    fn index(&self, i: usize) -> &Self::Output {
        return &self._i[i];
    }
}


#[derive(Debug)]
pub struct Sphere {
    pub r: f32,
    pub transform: Matrix,
    pub material: Material
}

impl Sphere {

    pub fn new(r:f32) -> Sphere {
        Sphere::new_with_transform(r, Matrix::identity())
    }

    pub fn new_with_transform(r:f32, m:Matrix) -> Sphere {
        Sphere {
            r:r,
            transform:m.inverse(),
            material:Default::default()
        }
    }

    pub fn new_with_transform_and_material(r:f32, m:Matrix, mat:Material) -> Sphere {
        Sphere {
            r:r,
            transform:m.inverse(),
            material:mat
        }
    }
}

impl Default for Sphere {

    fn default() -> Self {
        Sphere { r:1.0, transform:Matrix::identity(), material:Default::default() }
    }
}

impl Shape for Sphere {

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

    fn normal_at(&self, p: &Point) -> Vector {
        let object_point = p.transform(&self.transform);
        let object_normal = object_point.sub(Point::ZERO);
        let world_normal = object_normal.transform(&self.transform.transpose());
        return world_normal.normalize();
    }

    fn get_material(&self) -> &Material {
        &self.material
    }

    fn get_material_mut(&mut self) -> &mut Material {
        &mut self.material
    }
}

pub struct World {
    shapes: Vec<Box<dyn Shape>>,
    light: Light
}


impl World {

    pub fn new( light: Light) -> World {
        World { shapes: vec!(), light: light }
    }

    pub fn add_shape(& mut self, s: Box<dyn Shape>) -> & mut Self {
        self.shapes.push(s);
        self
    }

    pub fn get_shape(&self, i:usize) -> &dyn Shape {
        return self.shapes[i].borrow();
    }

    pub fn get_shape_mut<'a>(&'a mut self, i:usize) -> &'a mut (dyn Shape + 'static) {
        return self.shapes[i].borrow_mut();
    }

    pub fn intersect(&self, ray:&Ray) -> Intersections {
        let mut is = Intersections { _i:Vec::new() };
        for s in self.shapes.iter() {
            is.merge(s.intersect(ray));
        }
        return is;
    }

    pub fn shade_hit(&self, shape:&dyn Shape, comp:&CachedVectors) -> Color {
        return shape.get_material().lighting(&self.light, &comp.point, &comp.eyev, &comp.normal);  
    }

    pub fn color_at(&self, r:&Ray) -> Color {
        let is = self.intersect(r);
        return match is.get_hit() {
            Some(hit) => {
                self.shade_hit(hit.shape, &hit.compute_vectors(r))
            }
            None => { Color::BLACK }
        };
    }
}

impl Default for World {

    fn default() -> Self {
        let mut w = World { shapes:Vec::new(), light:Default::default() };
        let mut m = Material::default_with_color(Color::new(0.8, 1.0, 0.6));
        m.diffuse = 0.7;
        m.specular = 0.2;
        w.add_shape(Box::new(Sphere::new_with_transform_and_material(1.0, Matrix::identity(), m)));
        w.add_shape(Box::new(Sphere::new_with_transform(1.0, Matrix::identity().scaling(0.5, 0.5, 0.5))));
        return w;
    }
}


#[cfg(test)]
mod tests {

    use crate::vec::{Vector, Ray, Point};
    use super::*;
    use crate::*;
    
    #[test]
    fn test_intersect1() {
        let c = Sphere::new(1.);
        let ray = Ray::new(Point::new(0., 0., -5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);
        assert_eq!(v.len(), 2);
        
        assert!(floats_equal(v[0].t, 4.) || floats_equal(v[0].t, 6.));
        assert!(floats_equal(v[1].t, 4.) || floats_equal(v[1].t, 6.));

    }

    #[test]
    fn test_intersect2() {
        let c = Sphere::new(1.);
        let ray = Ray::new(Point::new(0., 1., -5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert!(floats_equal(v[0].t, 5.) && floats_equal(v[0].t, 5.));
        
    }

    #[test]
    fn test_intersect3() {
        let c = Sphere::new(1.);
        let ray = Ray::new(Point::new(0., 2., -5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert_eq!(v.len(), 0);
    }

    #[test]
    fn test_intersect4() {
        let c = Sphere::new(1.);
        let ray = Ray::new(Point::new(0., 0., 0.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert_eq!(v.len(), 2);
        
        assert!(floats_equal(v[0].t, -1.) || floats_equal(v[0].t, 1.));
        assert!(floats_equal(v[1].t, -1.) || floats_equal(v[1].t, 1.));
    }

    #[test]
    fn test_intersect5() {
        let c = Sphere::new(1.);
        let ray = Ray::new(Point::new(0., 0., 5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert_eq!(v.len(), 2);
        
        assert!(floats_equal(v[0].t, -4.) || floats_equal(v[0].t, -6.));
        assert!(floats_equal(v[1].t, -4.) || floats_equal(v[1].t, -6.));
    }

    #[test]
    fn test_hits1() {
        let c = Sphere::new(1.);
        let i1 = Intersection::new(1., &c);
        let i2 = Intersection::new(2., &c);
        let mut is = Intersections::new(i1);
        is.add(i2);

        assert_eq!(is.get_hit().unwrap().t, 1.);
    }

    #[test]
    fn test_hits2() {
        let c = Sphere::new(1.);
        let i1 = Intersection::new(-1., &c);
        let i2 = Intersection::new(-2., &c);
        let mut is = Intersections::new(i1);
        is.add(i2);

        assert!(is.get_hit().is_none());
    }

    #[test]
    fn test_shape_transform1() {
        let c = Sphere::new_with_transform(1., Matrix::identity().scaling(2., 2., 2.));
        let r = Ray::new(Point::new(0., 0., -5.),
                         Vector::new(0., 0., 1.));
        let is = c.intersect(&r);
        assert_eq!(2, is.len());
        
        assert!(floats_equal(is[0].t, 3.) || floats_equal(is[0].t, 7.));
        assert!(floats_equal(is[1].t, 3.) || floats_equal(is[1].t, 7.));
    }

    #[test]
    fn test_shape_transform2() {
        let c = Sphere::new_with_transform(1., Matrix::identity().translation(5., 0., 0.));
        let r = Ray::new(Point::new(0., 0., -5.),
                         Vector::new(0., 0., 1.));
        let is = c.intersect(&r);
        assert_eq!(0, is.len());
        
    }

    #[test]
    fn test_circle_normal1() {
        let c = Sphere::new(1.);
        let p = c.normal_at(&Point::new(1., 0., 0.));
        assert!(floats_equal(p.x, 1.));
        assert!(floats_equal(p.y, 0.));
        assert!(floats_equal(p.z, 0.));
    }

    #[test]
    fn test_circle_normal2() {
        let c = Sphere::new(1.);
        let p = c.normal_at(&Point::new(0., 1., 0.));
        assert!(floats_equal(p.x, 0.));
        assert!(floats_equal(p.y, 1.));
        assert!(floats_equal(p.z, 0.));
    }

    #[test]
    fn test_circle_normal3() {
        let c = Sphere::new(1.);
        let p = c.normal_at(&Point::new(0., 0., 1.));
        assert!(floats_equal(p.x, 0.));
        assert!(floats_equal(p.y, 0.));
        assert!(floats_equal(p.z, 1.));
    }

    #[test]
    fn test_circle_normal4() {
        let n = (3f32).sqrt()/3.0;
        
        let c = Sphere::new(1.);
        let p = c.normal_at(&Point::new(n, n, n));
        assert!(floats_equal(p.x, n));
        assert!(floats_equal(p.y, n));
        assert!(floats_equal(p.z, n));

    }
        
    #[test]
    fn test_circle_normal5() {
        let n = (3f32).sqrt()/3.0;
        
        let c = Sphere::new(1.);
        let p = c.normal_at(&Point::new(n, n, n));
        assert!(floats_equal(p.magnitude(), 1.))
    }

    #[test]
    fn test_translate_normal1() {
        let s = Sphere::new_with_transform(1., Matrix::identity().translation(0., 1., 0.));
        
        let v = s.normal_at(&Point::new(0., 1.70711, -0.70711));

        assert!(floats_equal(v.x, 0.));
        assert!(floats_equal(v.y, 0.70711));
        assert!(floats_equal(v.z, -0.70711));
                
    }

    #[test]
    fn test_translate_normal2() {
        // this is the same as scaling * rotation
        let s = Sphere::new_with_transform(1., Matrix::identity().
                                           rotation_z(std::f32::consts::PI/5.).
                                           scaling(1., 0.5, 1.));

        let n = (2f32).sqrt()/2.;
        
        let v = s.normal_at(&Point::new(0., n, -n));
        assert!(floats_equal(v.x, 0.));
        assert!(floats_equal(v.y, 0.97014));
        assert!(floats_equal(v.z, -0.24254));
        
    }

    #[test]
    fn test_world1() {
        let world:World = Default::default();
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let is = world.intersect(&ray);
        assert_eq!(is.len(), 4);
        let mut it = is.get_intersections_sorted();

        assert_eq!(it.next().unwrap().t, 4.);
        assert_eq!(it.next().unwrap().t, 4.5);
        assert_eq!(it.next().unwrap().t, 5.5);
        assert_eq!(it.next().unwrap().t, 6.);
        
    }

    #[test]
    fn test_world2() {
        let world:World = Default::default();
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let is = world.intersect(&ray);
        assert_eq!(is.len(), 4);
        
        assert_eq!((&is).into_iter().map(|a| a.t).collect::<Vec<f32>>(), vec![4., 4.5, 5.5, 6.]);
    }

    #[test]
    fn test_world3() {
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let s = Sphere::new(1.0);
        let is = s.intersect(&r);
        let ii:Vec<Intersection> = is.into_iter().filter(|a| a.t == 1.0).collect();
        let precomp = ii[0].compute_vectors(&r);
        assert_eq!(precomp.inside, true);
        assert_eq!(precomp.point, Point::new(0., 0., 1.));
        assert_eq!(precomp.eyev, Vector::new(0., 0., -1.));
    }

    #[test]
    fn test_world4() {
        let w:World = Default::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s:&dyn Shape = w.get_shape(0);
        let is = s.intersect(&r);
        let ii:Vec<Intersection> = is.into_iter().filter(|a| a.t == 4.0).collect();
        let precomp = ii[0].compute_vectors(&r);
        let c = w.shade_hit(s, &precomp);
        assert!(floats_equal(c.red, 0.38066));
        assert!(floats_equal(c.green, 0.47583));
        assert!(floats_equal(c.blue, 0.2855));
                
    }

    #[test]
    fn test_world5() {
        let mut w:World = Default::default();
        w.light = Light::new(Color::WHITE, Point::new(0., 0.25, 0.));
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let s:&dyn Shape = w.get_shape(1);
        let is = s.intersect(&r);
        let ii:Vec<Intersection> = is.into_iter().filter(|a| a.t == 0.5).collect();
        let precomp = ii[0].compute_vectors(&r);
        let c = w.shade_hit(s, &precomp);
        assert!(floats_equal(c.red, 0.90498));
        assert!(floats_equal(c.green, 0.90498));
        assert!(floats_equal(c.blue, 0.90498));
                
    }

    #[test]
    fn test_color_at1() {
        let w:World = Default::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 1., 0.));
        let c = w.color_at(&r);
        assert!(floats_equal(c.red, 0.));
        assert!(floats_equal(c.green, 0.));
        assert!(floats_equal(c.blue, 0.));
                
    }

    #[test]
    fn test_color_at2() {
        let w:World = Default::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let c = w.color_at(&r);
        assert!(floats_equal(c.red, 0.38066));
        assert!(floats_equal(c.green, 0.47583));
        assert!(floats_equal(c.blue, 0.2855));
                
    }

    #[test]
    fn test_color_at3() {
        let mut w:World = Default::default();
        let r = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        let color_compare;
        {
            let s = w.get_shape_mut(0);
            println!("outer color: {:?}", s.get_material().color);
            s.get_material_mut().ambient = 1.0;
            let s = w.get_shape_mut(1);
            println!("inner color: {:?}", s.get_material().color);
            s.get_material_mut().ambient = 1.0;
            color_compare = s.get_material().color;
        }
        let c = w.color_at(&r);
        assert!(floats_equal(c.red, color_compare.red));
        assert!(floats_equal(c.green, color_compare.green));
        assert!(floats_equal(c.blue, color_compare.blue));
                
    }

}

