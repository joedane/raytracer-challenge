

function initial_camera_position_circle (r, y)
   return {x=0, y=y, z=-r}
end

CIRCLE_SPAN=math.pi*2

function update_camera_position_circle (t, r, y)
   return {x=r*math.sin(t*CIRCLE_SPAN), y=y, z=-r*math.cos(t*CIRCLE_SPAN)}
end

function update_camera_position_linear (t) 
   local new_x = camera_start_position.x + (camera_end_position.x - camera_start_position.x) * t
   local new_y = camera_start_position.y + (camera_end_position.y - camera_start_position.y) * t
   local new_z = camera_start_position.z + (camera_end_position.z - camera_start_position.x) * t
   return {x=new_x, y=new_y, z=new_z}
end

function make_random_color()
   return { r=math.random(), g=math.random(), b=math.random() }
end

function make_random_scale() 
   return math.random()*2.0 + 0.50
end

function make_random_position()
   WORLD_WIDTH=10
   WORLD_HEIGHT=10
   return {
      x=math.random()*WORLD_WIDTH-(WORLD_WIDTH/2),
      y=math.random()*WORLD_HEIGHT,
      z=math.random()*WORLD_WIDTH-(WORLD_WIDTH/2),
   }
end

function make_random_shape ()
   local s = { type = "sphere", material=GLASS }
   s.color = make_random_color()
   s.scale = make_random_scale()
   s.position = make_random_position()
   return s
end
