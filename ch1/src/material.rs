
use crate::vec::{Point, Vector};
use crate::color::Color;
use crate::transform::Matrix;
use crate::shape::Shape;
use std::fmt;
use std::borrow::Borrow;


#[derive(Debug, Clone)]
pub struct Light {
    pub intensity: Color,
    pub position: Point
}

impl Light {

    pub fn new(i: Color, p: Point) -> Light {
        Light {
            intensity:i,
            position:p
        }
    }
}

impl Default for Light {

    fn default() -> Self {
        return Light::new(Color::WHITE, Point::new(-10., 10., -10.));
    }
}

pub trait Pattern : Sync {
    
    fn set_transform(&mut self, xf:Matrix);

    fn get_transform_inverse(&self) -> &Matrix;

    fn pattern_at(&self, point:&Point) -> Color;
    
    fn pattern_at_shape(&self, shape:&dyn Shape, world_point:&Point) -> Color {
        let object_point = world_point.transform(shape.get_transform_inverse());
        let pattern_point = object_point.transform(self.get_transform_inverse());
        return self.pattern_at(&pattern_point);
    }
}

pub struct TestPattern {
    xf_inv:Matrix,
}

impl TestPattern {

    pub fn new() -> Self {
        TestPattern { xf_inv:Matrix::identity() }
    }
}

impl Pattern for TestPattern {
    
    fn set_transform(&mut self, xf:Matrix) {
        self.xf_inv = xf.inverse();
    }
    
    fn get_transform_inverse(&self) -> &Matrix { &self.xf_inv }
    
    fn pattern_at(&self, point:&Point) -> Color {
        Color::new(point.x, point.y, point.z) // for testing purposes only
    }
}

pub struct StripePattern {
    xf_inv: Matrix,
    color_a: Color,
    color_b: Color,
}

impl StripePattern {
    
    pub fn new(a:Color, b:Color) -> Self {
        StripePattern {
            xf_inv: Matrix::identity(),
            color_a:a,
            color_b: b,
        }
    }
}

impl Pattern for StripePattern {

    fn set_transform(&mut self, xf:Matrix) {
        self.xf_inv = xf.inverse();
    }
    
    fn get_transform_inverse(&self) -> &Matrix { &self.xf_inv }
    
    fn pattern_at(&self, point:&Point) -> Color {
        if point.x.floor() % 2.0 == 0. {
            return self.color_a;
        } else {
            return self.color_b;
        }
    }
}

pub struct GradientPattern {
    xf_inv: Matrix,
    color_a: Color,
    color_b: Color,
    _diff: Color,
}

impl GradientPattern {
    
    pub fn new(a:Color, b:Color) -> Self {
        GradientPattern {
            xf_inv: Matrix::identity(),
            color_a:a,
            color_b: b,
            _diff: b.sub(a)
        }
    }
}

impl Pattern for GradientPattern {

    fn set_transform(&mut self, xf:Matrix) {
        self.xf_inv = xf.inverse();
    }
    
    fn get_transform_inverse(&self) -> &Matrix { &self.xf_inv }
    
    fn pattern_at(&self, point:&Point) -> Color {
        return self.color_a.add(self._diff.mul_f64(point.x - point.x.floor()));
    }
}

pub struct RingPattern {
    xf_inv: Matrix,
    color_a: Color,
    color_b: Color,
}

impl RingPattern {
    
    pub fn new(a:Color, b:Color) -> Self {
        RingPattern {
            xf_inv: Matrix::identity(),
            color_a:a,
            color_b: b,
        }
    }
}

impl Pattern for RingPattern {

    fn set_transform(&mut self, xf:Matrix) {
        self.xf_inv = xf.inverse();
    }
    
    fn get_transform_inverse(&self) -> &Matrix { &self.xf_inv }
    
    fn pattern_at(&self, point:&Point) -> Color {
        if (point.x.powi(2) + point.y.powi(2)).sqrt().floor() % 2.0 == 0.0 {
            return self.color_a;
        }
        else {
            return self.color_b;
        }
    }
}

pub struct CheckerPattern {
    xf_inv: Matrix,
    color_a: Color,
    color_b: Color,
}

impl CheckerPattern {
    
    pub fn new(a:Color, b:Color) -> Self {
        CheckerPattern {
            xf_inv: Matrix::identity(),
            color_a:a,
            color_b: b,
        }
    }
}

impl Pattern for CheckerPattern {

    fn set_transform(&mut self, xf:Matrix) {
        self.xf_inv = xf.inverse();
    }
    
    fn get_transform_inverse(&self) -> &Matrix { &self.xf_inv }
    
    fn pattern_at(&self, point:&Point) -> Color {
        if (point.x.floor() + point.y.floor() + point.z.floor()) % 2.0 == 0.0 {
            return self.color_a;
        }
        else {
            return self.color_b;
        }
    }
}

pub struct GridPattern {
    xf_inv: Matrix,
    color_base: Color,
    color_grid: Color,
}

impl GridPattern {
    
    pub fn new(base:Color, grid:Color) -> Self {
        GridPattern {
            xf_inv: Matrix::identity(),
            color_base:base,
            color_grid: grid,
        }
    }
}

impl Pattern for GridPattern {

    fn set_transform(&mut self, xf:Matrix) {
        self.xf_inv = xf.inverse();
    }
    
    fn get_transform_inverse(&self) -> &Matrix { &self.xf_inv }
    
    fn pattern_at(&self, point:&Point) -> Color {
        if (point.x - point.x.floor()).abs() < 0.01 ||
            (point.z - point.z.floor()).abs() < 0.01 {
                return self.color_grid;
            }
        else {
            return self.color_base;
        }
    }
}

pub struct Material {
    pub pattern: Option<Box<dyn Pattern>>,
    pub color: Option<Color>,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
    pub reflectiveness: f64,
    pub transparency:f64,
    pub refractive_index:f64,
}

impl fmt::Debug for Material {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Material").
            field("am", &self.ambient).
            field("di", &self.diffuse).
            field("sp", &self.specular).
            field("sh", &self.shininess).
            field("re", &self.reflectiveness).
            field("tr", &self.transparency).
            field("ri", &self.refractive_index).
            finish()
    }
}

impl Material {

    pub const DEFAULT:Material = Material {
        pattern: None,
        color: Some(Color::RED),
        ambient:0.1,
        diffuse:0.9,
        specular:0.9,
        shininess:200.0,
        reflectiveness:0.0,
        transparency:0.0,
        refractive_index:1.0,
    };
        

    pub fn solid(color:Color, ambient:f64, diffuse:f64, specular:f64, 
                 shininess:f64, reflectiveness:f64, transparency:f64, refractive_index:f64) -> Material {
        Material { 
            pattern:None, 
            color:Some(color),
            ambient, diffuse, specular, shininess, reflectiveness, transparency, refractive_index
        } 
    }

    pub fn solid_with_defaults(color:Color) -> Material {
        Material { color:Some(color), ..Material::DEFAULT } 
    }

    pub fn pattern_with_defaults(pattern:Box<dyn Pattern>) -> Material {
        Material { pattern:Some(pattern), ..Material::DEFAULT }
    }

    pub fn set_pattern(&mut self, pattern:Box<dyn Pattern>) -> &Self {
        self.pattern = Some(pattern);
        self
    }

    pub fn get_pattern(&self) -> &dyn Pattern {
        match &self.pattern {
            Some(p) => {
                p.borrow()
            }
            None => {
                panic!("no pattern set")
            }
        }
    }

    pub fn lighting(&self, light:&Light,
                    object:Option<&dyn Shape>,
                    point:&Point,
                    eye:&Vector,
                    normal:&Vector,
                    in_shadow: bool) -> Color {

        let effective_color = match &self.pattern {
            Some(p) => {
                p.pattern_at_shape(object.expect("must pass a shape when using a pattern"), point).mul(light.intensity)
            }
            None => {
                self.color.expect("must specify a color on materials when no pattern is used").mul(light.intensity)
            }
        };

        let lightv = light.position.sub(*point).normalize();
        let ambient = effective_color.mul_f64(self.ambient);

        if in_shadow {
            return ambient;
        }
        else {
            let light_dot_normal = lightv.dot(normal);
            let (diffuse, specular);
        
            if light_dot_normal < 0. {
                diffuse = Color::BLACK;
                specular = Color::BLACK;
            } else {
                diffuse = effective_color.mul_f64(self.diffuse * light_dot_normal);
                let reflectv = lightv.negate().reflect(normal);
                let reflect_dot_eye = reflectv.dot(eye);
                if reflect_dot_eye <= 0. {
                    specular = Color::BLACK;
                } else {
                    let factor = reflect_dot_eye.powf(self.shininess);
                    specular = light.intensity.mul_f64(self.specular * factor)
                }
            }
            return ambient.add(diffuse).add(specular);
        }
    }
}

impl Default for Material {

    fn default() -> Self {
        Material::solid_with_defaults(Color::WHITE)
    }
}

#[cfg(test)]
mod tests {

    use crate::vec::{Vector, Point};
    use crate::shape::{World, Sphere};
    use super::*;
    use crate::*;
    
    #[test]
    fn test_lighting1() {
        let m:Material = Default::default();
        let position = Point::ZERO;
        let eye = Vector::new(0., 0., -1.);
        let normal = Vector::new(0., 0., -1.);
        let light = Light::new(Color::new(1., 1., 1.), Point::new(0., 0., -10.));
        let result = m.lighting(&light, None, &position, &eye, &normal, false);
        assert!(floats_equal(result.red, 1.9));
        assert!(floats_equal(result.green, 1.9));
        assert!(floats_equal(result.blue, 1.9));
        
    }

    #[test]
    fn test_lighting2() {
        let m:Material = Default::default();
        let position = Point::ZERO;
        let n = (2f64).sqrt()/2.;
        let eye = Vector::new(0., n, -n);
        let normal = Vector::new(0., 0., -1.);
        let light = Light::new(Color::new(1., 1., 1.), Point::new(0., 0., -10.));
        let result = m.lighting(&light, None, &position, &eye, &normal, false);
        assert!(floats_equal(result.red, 1.));
        assert!(floats_equal(result.green, 1.));
        assert!(floats_equal(result.blue, 1.));
        
    }

    #[test]
    fn test_lighting3() {
        let m:Material = Default::default();
        let position = Point::ZERO;
        let eye = Vector::new(0., 0., -1.);
        let normal = Vector::new(0., 0., -1.);
        let light = Light::new(Color::new(1., 1., 1.), Point::new(0., 10., -10.));
        let result = m.lighting(&light, None, &position, &eye, &normal, false);
        assert!(floats_equal(result.red, 0.7364));
        assert!(floats_equal(result.green, 0.7364));
        assert!(floats_equal(result.blue, 0.7364));
    }

    #[test]
    fn test_lighting4() {
        let m:Material = Default::default();
        let position = Point::ZERO;
        let n = (2f64).sqrt()/2.;
        let eye = Vector::new(0., -n, -n);
        let normal = Vector::new(0., 0., -1.);
        let light = Light::new(Color::new(1., 1., 1.), Point::new(0., 10., -10.));
        let result = m.lighting(&light, None, &position, &eye, &normal, false);
        assert!(floats_equal(result.red, 1.6364));
        assert!(floats_equal(result.green, 1.6364));
        assert!(floats_equal(result.blue, 1.6364));
    }

    #[test]
    fn test_lighting5() {
        let m:Material = Default::default();
        let position = Point::ZERO;
        let eye = Vector::new(0., 0., -1.);
        let normal = Vector::new(0., 0., -1.);
        let light = Light::new(Color::new(1., 1., 1.), Point::new(0., 0., 10.));
        let result = m.lighting(&light, None, &position, &eye, &normal, false);
        println!("lighting5: {:?}", result);
        assert!(floats_equal(result.red, 0.1));
        assert!(floats_equal(result.green, 0.1));
        assert!(floats_equal(result.blue, 0.1));
    }

    #[test]
    fn test_shadow1() {
        let m:Material = Default::default();
        let eye = Vector::new(0., 0., -1.);
        let normal = Vector::new(0., 0., -1.);
        let light = Light::new(Color::WHITE, Point::new(0., 0., -10.));
        let color = m.lighting(&light, None, &Point::new(0., 0., 0.), &eye, &normal, true);
        assert_eq!(color, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn test_shadow2() {
        let w:World = Default::default();
        assert_eq!(w.is_shadowed(&Point::new(0., 10., 0.)), false);
    }

    #[test]
    fn test_shadow3() {
        let w:World = Default::default();
        assert_eq!(w.is_shadowed(&Point::new(10., -10., 10.)), true);
    }

    #[test]
    fn test_shadow4() {
        let w:World = Default::default();
        assert_eq!(w.is_shadowed(&Point::new(-20., 20., -20.)), false);
    }
 
    #[test]
    fn test_shadow5() {
        let w:World = Default::default();
        assert_eq!(w.is_shadowed(&Point::new(-2., 2., -2.)), false);
    }
    
    #[test]
    fn test_pattern1() {
        let mut material:Material = Default::default();
        let pattern = TestPattern::new();
        material.set_pattern(Box::new(pattern));
        let shape = Sphere::new_with_transform_and_material
        (Matrix::identity().scaling(2., 2., 2.), material);
        
        let c = shape.get_material().get_pattern().pattern_at_shape(&shape, &Point::new(2., 3., 4.));
        assert_eq!(c, Color::new(1., 1.5, 2.));
    }

    #[test]
    fn test_pattern2() {
        let mut material:Material = Default::default();
        let mut pattern = TestPattern::new();
        pattern.set_transform(Matrix::identity().scaling(2., 2., 2.));
        material.set_pattern(Box::new(pattern));
        let mut shape = Sphere::new();
        shape.material = material;
        let c = shape.get_material().get_pattern().pattern_at_shape(&shape, &Point::new(2., 3., 4.));
        assert_eq!(c, Color::new(1., 1.5, 2.));
    }

    #[test]
    fn test_pattern3() {
        let mut material:Material = Default::default();
        let mut pattern = TestPattern::new();
        pattern.set_transform(Matrix::identity().translation(0.5, 1., 1.5));
        material.set_pattern(Box::new(pattern));
        let shape = Sphere::new_with_transform_and_material
            ( Matrix::identity().scaling(2., 2., 2.), material);
        
        let c = shape.get_material().get_pattern().pattern_at_shape(&shape, &Point::new(2.5, 3., 3.5));
        assert_eq!(c, Color::new(0.75, 0.5, 0.25));
    }

    #[test]
    fn test_stripe_pattern1() {
        let mut material:Material = Default::default();
        let pattern = StripePattern::new(Color::WHITE, Color::BLACK);
        material.set_pattern(Box::new(pattern));
        let shape = Sphere::new_with_transform_and_material
            ( Matrix::identity().scaling(2., 2., 2.), material);
        
        let c = shape.get_material().get_pattern().pattern_at_shape(&shape, &Point::new(1.5, 0., 0.));
        assert_eq!(c, Color::WHITE);
    }

    #[test]
    fn test_stripe_pattern2() {
        let mut material:Material = Default::default();
        let mut pattern = StripePattern::new(Color::WHITE, Color::BLACK);
        pattern.set_transform(Matrix::identity().scaling(2., 2., 2.));
        material.set_pattern(Box::new(pattern));
        let mut shape = Sphere::new();
        shape.material = material;
        let c = shape.get_material().get_pattern().pattern_at_shape(&shape, &Point::new(1.5, 0., 0.));
        assert_eq!(c, Color::WHITE);
    }

    #[test]
    fn test_stripe_pattern3() {
        let mut material:Material = Default::default();
        let mut pattern = StripePattern::new(Color::WHITE, Color::BLACK);
        pattern.set_transform(Matrix::identity().translation(0.5, 0., 0.));
        material.set_pattern(Box::new(pattern));
        let shape = Sphere::new_with_transform_and_material
            (Matrix::identity().scaling(2., 2., 2.), material);
        
        let c = shape.get_material().get_pattern().pattern_at_shape(&shape, &Point::new(2.5, 0., 0.));
        assert_eq!(c, Color::WHITE);
    }

    #[test]
    fn test_gradient_pattern() {
        let p = GradientPattern::new(Color::WHITE, Color::BLACK);
        assert_eq!(p.pattern_at(&Point::new(0., 0., 0.)), Color::WHITE);
        assert_eq!(p.pattern_at(&Point::new(0.25, 0., 0.)), Color::new(0.75, 0.75, 0.75));
        assert_eq!(p.pattern_at(&Point::new(0.5, 0., 0.)), Color::new(0.5, 0.5, 0.5));
        assert_eq!(p.pattern_at(&Point::new(0.75, 0., 0.)), Color::new(0.25, 0.25, 0.25));
        
    }
}
