import uigen

ui.shader({
	'name':'alColorSpace',
	'description':'Converts a color to linear rec709 from a given color space',
	'output':'rgb',
	'maya_name':'alColorSpace',
	'maya_classification':'texture',
	'maya_id':'0x00116411',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_ColorSpace',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('input', 'rgb')
ui.parameter('sourceSpace', 'enum', 'Source space', enum_names=['sRGB', 'Cineon', 'LogC'])