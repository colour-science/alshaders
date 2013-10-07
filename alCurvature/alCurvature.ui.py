import uigen

ui.shader({
	'name':'alCurvature',
	'description':'Provides a measure of the surface curvature',
	'output':'rgb',
	'maya_name':'alCurvature',
	'maya_classification':'texture',
	'maya_id':'0x00116412',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_Curvature',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('samples', 'int')
ui.parameter('sampleOffset', 'float')
ui.parameter('sampleRadius', 'float')

uigen.remapControls(ui)

ui.parameter('color1', 'rgb')
ui.parameter('color1', 'rgb')