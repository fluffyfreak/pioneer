define_model('sphere', {
	info = {
 			lod_pixels = {2, 5, 10, 0},
			bounding_radius = 2,
			materials = {'body'},
			},

	static = function(lod)
	set_material('body', .5,.5,.5,1)
	use_material('body')
	sphere(2)
	end
})
