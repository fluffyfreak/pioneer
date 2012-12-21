define_model('halfSphere', {
	info = {
 			lod_pixels = {2, 5, 10, 0},
			bounding_radius = 0.5,
			materials = {'body'},
			},

	static = function(lod)
	set_material('body', .5,.5,.5,1)
	use_material('body')
	sphere_slice(6,6,0.5*math.pi,0.0, Matrix.scale(v(0.5,0.5,0.5))) 
	end
})
