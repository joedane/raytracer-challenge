
require ("functions")

WHITE = {r=1, g=1, b=1}
RED = {r=1, g=0, b=0}
GREEN = { r=0, g=1, b=0 }
BLUE = {r=0, g=0, b=1}
BLUISH = {r=0.537, g=0.831, b=0.914}

GLASS = { 
   ambient=0.1, diffuse=0.5, 
   specular=0.3, shininess=100.0,
   reflectiveness=0.5,transparency=1.0,
   refractive_index=1.1
}

SOLID = {
   ambient=0.1, diffuse=0.7, 
   specular=0.0, shininess=100.0,
   reflectiveness=0.1,transparency=0,
}

CHECKS = { type = "checks", scale=4.0 }
GRID = { type = "grid", scale = 1.0 }

world = {
   lights = {
      { color = WHITE, position = {x=0, y=10, z=-10} },
      { color = RED, position = {x=5, y=10, z=10} }
   },
   shapes = {
--      { type = "plane", color = {r=179/255, g=230/255, b=255/255}, y=-3.0 },
      { type = "sphere", material = GLASS, color={0,0,0}, scale = 1.5, position = {x=-3.5, y=0, z=0} }, 
      { type = "cube", material=SOLID, color=RED, rotate_y=45 },
      { type = "sphere", material = SOLID, color = GREEN, position = {x=1.5, y=0, z=3.5} },
      --      { type = "plane_xz", pattern = CHECKS, y=-3.0 }
   }
}

camera_start_position = {x=0, y=5, z=-20}
camera_end_position = {x=-10, y=1, z=-10}

camera = {
   screenwidth = 600,
   screenheight = 400,
   position = initial_camera_position_circle(10, 3),
   lookat = {x=0, y=0, z=0},
   up = {x=0, y=1, z=0},
   fov = math.pi/3,
   samples = 50,
}


function populate_world (w, num_shapes) 
   for i=1, num_shapes, 1 do
      local s = make_random_shape()
      print(string.format("adding shape at (%f, %f, %f)", s.position.x, s.position.y, s.position.z))
      table.insert(w.shapes, s)
   end
end

function run1() 
   math.randomseed(13)
--   populate_world(world, 10)
   Render(world, camera, "test.jpg")
end

function run ()
   math.randomseed(13)
   populate_world(world, 10)
   local out = StartAnimation("anim2.gif")
   local n_frames = 30
   for i=1, n_frames, 1 do
      out:AddFrame(world, camera)
      print (string.format("Frame %d complete", i)) 
      camera.position = update_camera_position_circle(i/n_frames, 10, 5)
      print (string.format("moved to (x,z) (%f, %f)", camera.position.x, camera.position.z))
   end
   out:Finish()
end

run()
