import uigen

ui.shader({
	'name':'alInputVector',
	'description':'Grab useful vectors from the shader state and use them in a network',
	'output':'vector',
	'maya_name':'alInputVector',
	'maya_classification':'texture',
	'maya_id':'0x0011640E',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_InputVector',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('input', 'enum', 'Input', enum_names=[
	"P",
	"Po",
	"N",
	"Nf",
	"Ng",
	"Ngf",
	"Ns",
	"dPdu",
	"dPdv",
	"Ld",
	"Rd",
	"uv",
	"User",
	"Custom"
])

ui.parameter('userName', 'string', 'User name')
ui.parameter('vector', 'vector', 'Custom')
ui.parameter('type', 'enum', 'Type', enum_names=[
	"Point",
	"Vector"
])

with uigen.group(ui, 'Transform'):
	ui.parameter('matrix', 'matrix', 'Matrix')
	ui.parameter('coordinates', 'enum', 'Coordinates', enum_names=[
		"cartesian",
		"spherical",
		"normalized spherical"
	])