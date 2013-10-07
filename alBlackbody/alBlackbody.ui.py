import uigen

ui.shader({
	'name':'alBlackbody',
	'description':'Generates a rec709 color from the blackbody spectrum for the given temperature',
	'output':'rgb',
	'maya_name':'alBlackbody',
	'maya_classification':'texture',
	'maya_id':'0x00116404',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_Blackbody',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('temperature', 'float', 1000.0, mn=273, mx=16000)
ui.parameter('strength', 'float', 1.0)

with uigen.group(ui, 'Advanced'):
	ui.parameter('physicalIntensity', 'float', 1.0, 'Physical intensity')
	ui.parameter('physicalExposure', 'float', -16.0, 'Physical exposure')