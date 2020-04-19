
PI = 3.14159265358979323846
WHITE = {r=1, g=1, b=1}
RED = {r=1, g=0, b=0}
GREEN = { r=0, g=1, b=0 }

GLASS = { 
   ambient=0.0, diffuse=0.9, 
   specular=0.3, shininess=100.0,
   reflectiveness=1.0,transparency=1.0,
   refractive_index=1.6
}

CHECKS = { type = "checks", scale=4.0 }
GRID = { type = "grid", scale = 1.0 }

world = {
   lights = {
      { color = WHITE, position = {x=0, y=10, z=-10} },
      { color = RED, position = {x=5, y=10, z=10} }
   },
   shapes = {
      { type = "sphere", material = GLASS, color={0.1, 0.1, 0.1}, scale = 1.0, position = {x=0, y=0, z=0} }, 
      { type = "sphere", material = GLASS, color = GREEN, position = {x=1.5, y=0, z=3.5} },
      --            { type = "plane_xz", color = {r=179/255, g=230/255, b=255/255}, y=-3.0 }
      { type = "plane_xz", pattern = CHECKS, y=-3.0 }
   }
}

camera_start_position = {x=0, y=1, z=-10}
camera_end_position = {x=-10, y=1, z=-10}

camera = {
   screenwidth = 600,
   screenheight = 400,
   position = camera_start_position,
   lookat = {x=0, y=0, z=0},
   up = {x=0, y=1, z=0},
   fov = PI/3,
   samples = 1,
}

function update_camera_position (t) 
   local new_x = camera_start_position.x + (camera_end_position.x - camera_start_position.x) * t
   local new_y = camera_start_position.y + (camera_end_position.y - camera_start_position.y) * t
   local new_z = camera_start_position.z + (camera_end_position.z - camera_start_position.x) * t
   return {x=new_x, y=new_y, z=new_z}
end

function run ()

   local out = StartAnimation("anim2.gif")
   local n_frames = 1
   for i=1, n_frames, 1 do
      out:AddFrame(world, camera)
      print (string.format("Frame %d complete", i)) 
      camera.position = update_camera_position(i/n_frames)
   end
   out:Finish()
end

run()
