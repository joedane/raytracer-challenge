
use crate::vec::{Point, Vector};
use crate::color::Color;


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

pub trait Pattern {
    
    fn set_transform(&mut self, xf:Matrix);

    fn get_trasnform(&self) -> &Matrix;

    fn pattern_at(&self, shape:&dyn Shape, point:&Point) {
        
    }
}


#[derive(Debug, Clone)]
pub struct Material {
    pub color: Color,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64
}

impl Material {

    pub const DEFAULT:Material = Material {
        color:Color::RED, 
        ambient:0.1,
        diffuse:0.9,
        specular:0.9,
        shininess:200.0
    };
        

    pub fn new(color:Color, ambient:f64, diffuse:f64, specular:f64, shininess:f64) -> Material {
        Material { color, ambient, diffuse, specular, shininess }
    }

    pub fn default_with_color(c:Color) -> Material {
        Material { color:c, ..Material::DEFAULT }
    }

    pub fn lighting(&self, light:&Light, 
                    point:&Point,
                    eye:&Vector,
                    normal:&Vector,
                    in_shadow: bool) -> Color {

        let effective_color = self.color.mul(light.intensity);
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
        Material::new(Color::WHITE, 0.1, 0.9, 0.9, 200.0)
    }
}

#[cfg(test)]
mod tests {

    use crate::vec::{Vector, Point};
    use crate::shape::World;
    use super::*;
    use crate::*;
    
    #[test]
    fn test_lighting1() {
        let m:Material = Default::default();
        let position = Point::ZERO;
        let eye = Vector::new(0., 0., -1.);
        let normal = Vector::new(0., 0., -1.);
        let light = Light::new(Color::new(1., 1., 1.), Point::new(0., 0., -10.));
        let result = m.lighting(&light, &position, &eye, &normal, false);
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
        let result = m.lighting(&light, &position, &eye, &normal, false);
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
        let result = m.lighting(&light, &position, &eye, &normal, false);
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
        let result = m.lighting(&light, &position, &eye, &normal, false);
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
        let result = m.lighting(&light, &position, &eye, &normal, false);
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
        let color = m.lighting(&light, &Point::new(0., 0., 0.), &eye, &normal, true);
        assert_eq!(color, Color::new(0.1, 0.1, 0.1));
    }

    fn test_shadow2() {
        let w:World = Default::default();
        assert_eq!(w.is_shadowed(&Point::new(0., 10., 0.)), false);
    }

    fn test_shadow3() {
        let w:World = Default::default();
        assert_eq!(w.is_shadowed(&Point::new(10., -10., 10.)), true);
    }
 
    fn test_shadow4() {
        let w:World = Default::default();
        assert_eq!(w.is_shadowed(&Point::new(-20., 20., -20.)), false);
    }
 
    fn test_shadow5() {
        let w:World = Default::default();
        assert_eq!(w.is_shadowed(&Point::new(-2., 2., -2.)), true);
    }
  
}
