import uigen

ui.shader({
	'name':'alRemapColor',
	'description':'Controls to adjust a color value',
	'output':'rgb',
	'maya_name':'alRemapColor',
	'maya_classification':'texture',
	'maya_id':'0x0011640D',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_RemapColor',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('input', 'rgb', 'Input')
ui.parameter('gamma', 'float', 'Gamma')
ui.parameter('saturation', 'float', 'Saturation')
ui.parameter('hueOffset', 'float', 'Hue offset')

with uigen.group(ui, 'Contrast'):
	ui.parameter('contrast', 'float', 'Contrast')
	ui.parameter('contrastPivot', 'float', 'Pivot')

ui.parameter('gain', 'float', 'Gain')
ui.parameter('exposure', 'float', 'Exposure')
ui.parameter('mask', 'float', 'Mask')