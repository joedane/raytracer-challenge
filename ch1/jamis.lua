
FLOOR = {
   type = "checks"
}
WHITE = {r=1, g=1, b=1}

WALL_MATERIAL = {
   pattern = {
      type = "stripes",
      color_a = {r=0.45, g=0.45, b=0.45},
      color_b = {r=0.55, g=0.55, b=0.55},
      scale = 0.25,
      rotate_y = 1.5708,   
   },
   ambient = 0,
   diffuse = 0.4,
   specular = 0,
   reflective = 0.3
}

camera = {
   screenwidth = 600,
   screenheight = 400,
   position = {x=-2.6, y=1.5, z=-3.9},
   lookat = {x=-0.6, y=1, z=-0.8},
   up = {x=0, y=1, z=0},
   fov = 1.152,
   samples = 50,
}

world = {
   lights = {
      { color = WHITE, position={x=-4.9, y=4.9, z=-1} },
   },
   shapes = {
      {  type = "plane", 
         rotate_y = 0.31415, 
         material = { 
            pattern = {  type = "checks", color_a = {r=0.35, g=0.35, b=0.35}, color_b = {r=0.65, g=0.65, b=0.65}},
            specular = 0,
            reflectiveness = 0.4         
         }
      },
      { type = "sphere", scale=0.4, position={x=4.6, y=0.4, z=1}, material = { color={r=0.8, g=0.5, b=0.3}, shininess=50 }},
      { type = "sphere", scale=0.3, position={x=4.7, y=0.3, z=0.4}, material = { color={r=0.9, g=0.4, b=0.5}, shininess=50 }},
      { type = "sphere", scale=0.5, position={x=-1, y=0.5, z=4.5}, material = { color={r=0.4, g=0.9, b=0.6}, shininess=50 }},
      { type = "sphere", scale=0.3, position={x=-1.7, y=0.3, z=4.7}, material = { color={r=0.4, g=0.6, b=0.9}, shininess=50 }}, 
      { type = "sphere", scale=1, position={x=-0.6, y=1, z=0.6}, material = { color={r=1, g=0.3, b=0.2}, shininess=5, specular=0.4 }},
      { type = "sphere", scale=0.7, position={x=0.6, y=0.7, z=-0.6}, material = { color={r=0, g=0, b=0.2}, ambient=0, diffuse=0.4, specular=0.9, shininess=300, reflectiveness=0.9, transparency=0.9, refractive_index=1.5}}, 
      { type = "sphere", scale=0.5, position={x=-0.7, y=0.5, z=-0.8}, material = { color={r=0, g=0.2, b=0}, ambient=0, diffuse=0.4, specular=0.9, shininess=300, reflectiveness=0.9, transparency=0.9, refractive_index=1.5}}, 

   }
}


Render(world, camera, "jamis.jpg")
