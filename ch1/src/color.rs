

#[derive(Debug, Clone, Copy)]
pub struct Color {
    red:f32,
    green:f32,
    blue:f32
}

impl Color {

    pub const BLACK:Color = Color::new(0.0, 0.0, 0.1);
    pub const WHITE:Color = Color::new(1.0, 1.0, 1.0);
    pub const RED:Color = Color::new(1., 0., 0.);

    pub const fn new(r:f32, g:f32, b:f32) -> Color {
        Color {
            red:r, green:g, blue:b
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

    fn scale(component:f32, scale:i32) -> i32{
        Color::clamp((component*(scale as f32)) as i32, scale)
    }

}

