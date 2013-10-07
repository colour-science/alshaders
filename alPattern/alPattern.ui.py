import uigen

ui.shader({
	'name':'alPattern',
	'description':'Simple pattern generator',
	'output':'rgb',
	'maya_name':'alPattern',
	'maya_classification':'texture',
	'maya_id':'0x0011640A',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_Pattern',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('space', 'enum', 'Space', enum_names=[
	"world",
	"object",
	"Pref",
	"UV"
])

ui.parameter('axis', 'enum', 'Axis', enum_names=[
	"X",
	"Y",
	"Z"
])

ui.parameter('shape', 'enum', 'Shape', enum_names=[
	"sine",
	"square",
	"saw"
])

ui.parameter('frequency', 'float', 'Frequency', connectible=False)
ui.parameter('offset', 'float', 'Offset', connectible=False)

uigen.remapControls(ui)

ui.parameter('color1', 'rgb', 'Color 1')
ui.parameter('color2', 'rgb', 'Color 2')
