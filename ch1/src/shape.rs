

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

    fn intersect(&self, ray: &Ray) -> Intersections {
        let ray = ray.transform(self.get_transform_inverse());
        return self.intersect_local(&ray);
    }

    fn intersect_local(&self, ray:&Ray) -> Intersections;
    
    fn get_transform_inverse(&self) -> &Matrix;
    
    fn normal_at(&self, p: &Point) -> Vector {
        let xf = &self.get_transform_inverse();
        let local_point = p.transform(xf);
        let local_normal = self.normal_at_local(&local_point);
        let world_normal = local_normal.transform(&xf.transpose());
        return world_normal.normalize();
    }

    fn normal_at_local(&self, p: &Point) -> Vector;

    fn get_material(&self) -> &Material;

    fn get_material_mut(&mut self) -> &mut Material;

    fn box_clone(&self) -> Box<dyn Shape>;
}

#[derive(Debug)]
pub struct CachedVectors {
    t:f64, 
    point:Point,
    over_point:Point,
    eyev:Vector,
    normal:Vector,
    inside:bool
}

impl CachedVectors {
    
    pub fn new(t:f64, point:Point, eyev:Vector, normal:Vector) -> Self{
        let inside = normal.dot(&eyev) < 0.0;
        let normal = if inside { normal.negate() } else { normal };
        let over_point = point.add_vec(normal.mul(Vector::EPSILON));
        
        return CachedVectors { 
            t:t, 
            point:point, 
            over_point:over_point,
            eyev:eyev, 
            normal: normal,
            inside:inside
        };
    }
}

#[derive(Debug)]
pub struct Intersection<'a> {
    pub t:f64,
    pub shape:&'a dyn Shape
}

impl<'a> Intersection<'a> {

    fn new(t:f64, shape:&'a dyn Shape) -> Intersection<'a> {
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

    pub fn new_empty() -> Intersections<'a> {
        Intersections {
            _i:Vec::new()
        }
    }
    
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


#[derive(Debug, Clone)]
pub struct Sphere {
    pub r: f64,
    pub transform: Matrix,
    pub material: Material,
    _force_ctor: (),  // see https://words.steveklabnik.com/structure-literals-vs-constructors-in-rust
}

impl Sphere {

    pub fn new() -> Sphere {
        Sphere::new_with_transform(1.0, Matrix::identity())
    }

    pub fn new_with_transform(r:f64, m:Matrix) -> Sphere {
        Sphere {
            r:r,
            transform:m.inverse(),
            material:Default::default(),
            _force_ctor: ()
        }
    }

    pub fn new_with_transform_and_material(r:f64, m:Matrix, mat:Material) -> Sphere {
        Sphere {
            r:r,
            transform:m.inverse(),
            material:mat,
            _force_ctor: ()
        }
    }
}

impl Default for Sphere {

    fn default() -> Self {
        Sphere { r:1.0, transform:Matrix::identity(), material:Default::default(), _force_ctor: () }
    }
}

impl Shape for Sphere {

    fn get_transform_inverse(&self) -> &Matrix {
        return &self.transform;
    }

    /**
     * the Ray has already been transformed to local coordinates.
     */
    fn intersect_local(&self, ray: &Ray) -> Intersections {

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

    fn normal_at_local(&self, p: &Point) -> Vector {
        return p.sub(Point::ZERO);
    }

    fn get_material(&self) -> &Material {
        &self.material
    }

    fn get_material_mut(&mut self) -> &mut Material {
        &mut self.material
    }

    fn box_clone(&self) -> Box<dyn Shape> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Shape> {

    fn clone(&self) -> Box<dyn Shape> {
        self.box_clone()
    }
}

#[derive(Clone, Debug)]
pub struct Plane {
    transform: Matrix,
    material: Material, 
    local_normal: Vector
}

impl Plane {

    pub fn new() -> Plane {
        Plane {
            transform: Matrix::identity(),
            material: Default::default(),
            local_normal: Vector::new(0., 1., 0.)
        }
    }
}

impl Shape for Plane {


    fn intersect_local(&self, ray:&Ray) -> Intersections {
        if ray.direction.y.abs() < Vector::EPSILON {
            return Intersections::new_empty();
        }
        else {
            let t = -ray.origin.y / ray.direction.y;
            let i = Intersection::new(t, self as &dyn Shape);
            return Intersections::new(i);
        }
    }
    
    fn get_transform_inverse(&self) -> &Matrix {
        &self.transform
    }
    
    fn normal_at_local(&self, _p: &Point) -> Vector {
        return self.local_normal;
    }

    fn get_material(&self) -> &Material {
        &self.material
    }

    fn get_material_mut(&mut self) -> &mut Material {
        &mut self.material
    }

    fn box_clone(&self) -> Box<dyn Shape> {
        Box::new(self.clone())
    }

}


#[derive(Clone)]
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

    pub fn shade_hit(&self, material:&Material, comp:&CachedVectors) -> Color {
        return material.lighting
            (&self.light, &comp.over_point, &comp.eyev, &comp.normal, self.is_shadowed(&comp.over_point));  
    }

    pub fn color_at(&self, r:&Ray) -> Color {
        let is = self.intersect(r);
        return match is.get_hit() {
            Some(hit) => {
                self.shade_hit(hit.shape.get_material(), &hit.compute_vectors(r))
            }
            None => { Color::BLACK }
        };
    }

    pub fn is_shadowed(&self, p:&Point) -> bool {
        let v = self.light.position.sub(*p);
        let distance = v.magnitude();
        let direction = v.normalize();
        let r = Ray::new(*p, direction);
        let is = self.intersect(&r);

        match is.get_hit() {
            Some(h) => { h.t < distance }
            None => { false }
        }
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
        let c = Sphere::new();
        let ray = Ray::new(Point::new(0., 0., -5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);
        assert_eq!(v.len(), 2);
        
        assert!(floats_equal(v[0].t, 4.) || floats_equal(v[0].t, 6.));
        assert!(floats_equal(v[1].t, 4.) || floats_equal(v[1].t, 6.));

    }

    #[test]
    fn test_intersect2() {
        let c = Sphere::new();
        let ray = Ray::new(Point::new(0., 1., -5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert!(floats_equal(v[0].t, 5.) && floats_equal(v[0].t, 5.));
        
    }

    #[test]
    fn test_intersect3() {
        let c = Sphere::new();
        let ray = Ray::new(Point::new(0., 2., -5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert_eq!(v.len(), 0);
    }

    #[test]
    fn test_intersect4() {
        let c = Sphere::new();
        let ray = Ray::new(Point::new(0., 0., 0.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert_eq!(v.len(), 2);
        
        assert!(floats_equal(v[0].t, -1.) || floats_equal(v[0].t, 1.));
        assert!(floats_equal(v[1].t, -1.) || floats_equal(v[1].t, 1.));
    }

    #[test]
    fn test_intersect5() {
        let c = Sphere::new();
        let ray = Ray::new(Point::new(0., 0., 5.), 
                           Vector::new(0., 0., 1.));
        let v = c.intersect(&ray);

        assert_eq!(v.len(), 2);
        
        assert!(floats_equal(v[0].t, -4.) || floats_equal(v[0].t, -6.));
        assert!(floats_equal(v[1].t, -4.) || floats_equal(v[1].t, -6.));
    }

    #[test]
    fn test_intersect_jd1() {
        let c = Sphere::new_with_transform(1.0, Matrix::identity().scaling(1.0000, 1.0000, 1.0000));

        let n = (std::f64::consts::PI*3.0/4.0).sin();
        
        let ray = Ray::new(Point::new(0., 0., -2.0*n), 
                           Vector::new(0., 1., 1.).normalize());
        let v = c.intersect(&ray);

        println!("{:?}", v.len());
        if v.len() == 2 {
            println!("{:?}", v[0]);
            println!("{:?}", v[1]);
        }
    }
    
    #[test]
    fn test_hits1() {
        let c = Sphere::new();
        let i1 = Intersection::new(1., &c);
        let i2 = Intersection::new(2., &c);
        let mut is = Intersections::new(i1);
        is.add(i2);

        assert_eq!(is.get_hit().unwrap().t, 1.);
    }

    #[test]
    fn test_hits2() {
        let c = Sphere::new();
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
        let c = Sphere::new();
        let p = c.normal_at(&Point::new(1., 0., 0.));
        assert!(floats_equal(p.x, 1.));
        assert!(floats_equal(p.y, 0.));
        assert!(floats_equal(p.z, 0.));
    }

    #[test]
    fn test_circle_normal2() {
        let c = Sphere::new();
        let p = c.normal_at(&Point::new(0., 1., 0.));
        assert!(floats_equal(p.x, 0.));
        assert!(floats_equal(p.y, 1.));
        assert!(floats_equal(p.z, 0.));
    }

    #[test]
    fn test_circle_normal3() {
        let c = Sphere::new();
        let p = c.normal_at(&Point::new(0., 0., 1.));
        assert!(floats_equal(p.x, 0.));
        assert!(floats_equal(p.y, 0.));
        assert!(floats_equal(p.z, 1.));
    }

    #[test]
    fn test_circle_normal4() {
        let n = (3f64).sqrt()/3.0;
        
        let c = Sphere::new();
        let p = c.normal_at(&Point::new(n, n, n));
        assert!(floats_equal(p.x, n));
        assert!(floats_equal(p.y, n));
        assert!(floats_equal(p.z, n));

    }
        
    #[test]
    fn test_circle_normal5() {
        let n = (3f64).sqrt()/3.0;
        
        let c = Sphere::new();
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
                                           rotation_z(std::f64::consts::PI/5.).
                                           scaling(1., 0.5, 1.));

        let n = (2f64).sqrt()/2.;
        
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
        
        assert_eq!((&is).into_iter().map(|a| a.t).collect::<Vec<f64>>(), vec![4., 4.5, 5.5, 6.]);
    }

    #[test]
    fn test_world3() {
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
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
        let c = w.shade_hit(s.get_material(), &precomp);
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
        let c = w.shade_hit(s.get_material(), &precomp);

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
            s.get_material_mut().ambient = 1.0;
            let s = w.get_shape_mut(1);
            s.get_material_mut().ambient = 1.0;
            color_compare = s.get_material().color;
        }
        let c = w.color_at(&r);
        assert!(floats_equal(c.red, color_compare.red));
        assert!(floats_equal(c.green, color_compare.green));
        assert!(floats_equal(c.blue, color_compare.blue));
                
    }

    #[test]
    fn test_color_at_jd1() {
        let mut w = World::new(Default::default());
        let mut c = Sphere::new_with_transform(1.0, Matrix::identity().scaling(1.00001, 1.00001, 1.00001));
        c.get_material_mut().color = Color::RED;
        w.add_shape(Box::new(c));
        
        let n = (std::f64::consts::PI*3.0/4.0).sin();
        
        let ray = Ray::new(Point::new(0., 0., -2.0*n), 
                           Vector::new(0., 1., 1.).normalize());

        let c = w.color_at(&ray);

        assert_ne!(c, Color::BLACK);

    }

    #[test]
    fn test_shadow1() {
        let mut w = World::new(Light::new(Color::WHITE, Point::new(0., 0., -10.)));
        w.add_shape(Box::new(Sphere::new()));
        w.add_shape(Box::new(Sphere::new_with_transform(1., Matrix::identity().translation(0., 0., 10.))));
        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let is = w.intersect(&r); 
        let ii:Vec<Intersection> = is.into_iter().filter(|a| a.t == 4.0).collect();       
        assert_eq!(ii.len(), 1);
        let vectors = ii[0].compute_vectors(&r);
        let c = w.shade_hit(ii[0].shape.get_material(), &vectors);
        assert_eq!(c, Color::new(0.1, 0.1, 0.1));        
    }

    #[test]
    fn test_shadow2() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let shape = Sphere::new_with_transform(1.0, Matrix::identity().translation(0., 0., 1.));
        let is = shape.intersect(&r);
        let ii:Vec<Intersection> = is.into_iter().filter(|a| a.t == 5.0).collect();       
        assert_eq!(ii.len(), 1);
        let vectors = ii[0].compute_vectors(&r);
        assert!(vectors.over_point.z < -Vector::EPSILON/2.0);
        assert!(vectors.point.z > vectors.over_point.z);
    }

    #[test]
    fn test_plane1() {
        let p = Plane::new();
        let (n1, n2, n3) = (p.normal_at_local(&Point::new(0., 0., 0.)),
                            p.normal_at_local(&Point::new(10., 0., -10.)),
                            p.normal_at_local(&Point::new(-5., 0., 150.)));
        assert_eq!(n1, Vector::new(0., 1., 0.));
        assert_eq!(n2, Vector::new(0., 1., 0.));
        assert_eq!(n3, Vector::new(0., 1., 0.));
    }

    #[test]
    fn test_plane2() {
        let p = Plane::new();
        let r = Ray::new(Point::new(0., 10., 0.), Vector::new(0., 0., 1.));
        let is = p.intersect(&r);
        assert_eq!(is.len(), 0);
    }

    #[test]
    fn test_plane3() {
        let p = Plane::new();
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let is = p.intersect(&r);
        assert_eq!(is.len(), 0);
    }

    #[test]
    fn test_plane4() {
        let p = Plane::new();
        let r = Ray::new(Point::new(0., 1., 0.), Vector::new(0., -1., 0.));
        let is = p.intersect(&r);
        assert_eq!(is.len(), 1);
        assert_eq!(is.into_iter().next().unwrap().t, 1.0)
    }


    #[test]
    fn test_plane5() {
        let p = Plane::new();
        let r = Ray::new(Point::new(0., -1., 0.), Vector::new(0., 1., 0.));
        let is = p.intersect(&r);
        assert_eq!(is.len(), 1);
        assert_eq!(is.into_iter().next().unwrap().t, 1.0)
            
    }

        
}

