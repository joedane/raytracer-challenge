
use std::i32;
use crate::floats_equal;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Color {
    pub red:f64,
    pub green:f64,
    pub blue:f64
}

impl Color {

    pub const BLACK:Color = Color::new(0.0, 0.0, 0.0);
    pub const WHITE:Color = Color::new(1.0, 1.0, 1.0);
    pub const RED:Color = Color::new(1., 0., 0.);
    pub const GREEN:Color = Color::new(0., 1., 0.);
    pub const BLUE:Color = Color::new(0., 0., 1.);

    pub const fn new(r:f64, g:f64, b:f64) -> Color {
        Color {
            red:r, green:g, blue:b
        }
    }

    fn color_component_from_spec(s:&str) -> Result<f64, String> {
        match i32::from_str_radix(s, 16) {
            Ok(v) => { Result::Ok(v as f64 / 255.0) }
            Err(e) => { Err(e.to_string()) }
        }
    }
    
    pub fn color_from_spec(s:&str) -> Result<Color, &'static str> {
        if s.len() != 7 || !s.starts_with('#') {
            Err("invalid color representaion")
        } else {
            Ok(Color::new(Color::color_component_from_spec(&s[1..3]).unwrap(),
                          Color::color_component_from_spec(&s[3..5]).unwrap(),
                          Color::color_component_from_spec(&s[5..7]).unwrap()))
        }
    }
    
    pub fn red_scaled(&self, scale:i32) -> i32 {
        Color::scale(self.red, scale)
    }

    pub fn green_scaled(&self, scale:i32) -> i32 {
        Color::scale(self.green, scale)
    }

    pub fn blue_scaled(&self, scale:i32) -> i32 {
        Color::scale(self.blue, scale)
    }

    pub fn add(&self, other:Color) -> Color {
        Color {
            red:self.red + other.red,
            green:self.green + other.green,
            blue:self.blue + other.blue
        }
    }
    
    pub fn sub(&self, other:Color) -> Color {
        Color {
            red:self.red - other.red,
            green:self.green - other.green,
            blue:self.blue - other.blue
        }
    }

    pub fn mul(&self, other:Color) -> Color {
        Color {
            red:self.red * other.red,
            green:self.green * other.green,
            blue:self.blue * other.blue
        }
    }
    
    // can this be combined with the above as a generic function?
    pub fn mul_f64(&self, m:f64) -> Color {
        Color {
            red:self.red * m,
            green:self.green * m,
            blue:self.blue * m
        }
    }

    fn clamp(i:i32, scale:i32) -> i32 {
        if i < 0 {
            0
        }
        else if i > scale {
            scale
        }
        else {
            i
        }
    }

    fn scale(component:f64, scale:i32) -> i32{
        Color::clamp((component*(scale as f64)) as i32, scale)
    }

    pub fn approximately_equal(&self, other:&Color) -> bool {
        floats_equal(self.red, other.red) &&
            floats_equal(self.green, other.green) &&
            floats_equal(self.blue, other.blue)
    }
}

#[cfg(test)]
mod tests {

    use super::Color;

    #[test]
    fn test_color1() {
        assert_eq!(Color::WHITE, Color::color_from_spec("#ffffff").unwrap());
        assert_eq!(Color::RED, Color::color_from_spec("#ff0000").unwrap());
    }
    
}
