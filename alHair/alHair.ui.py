import uigen

ui.shader({
	'name':'alHair',
	'description':'Dual-scattering, multiple-importance-sampled hair shader',
	'output':'rgb',
	'maya_name':'alHair',
	'maya_classification':'shader/surface',
	'maya_id':'0x00116403',
	'maya_swatch':True,
	'maya_matte':True,
	'maya_bump':False,
	'soft_name':'ALS_Hair',
	'soft_classification':'material',
	'soft_version':1
})

with uigen.group(ui, 'Fibre properties', False):
	ui.parameter('hairColor', 'rgb', 'Color')
	ui.parameter('specularWidth', 'float', 'Highlight width', connectible=False)
	ui.parameter('specularShift', 'float', 'Highlight shift', connectible=False)
	ui.parameter('opacity', 'rgb')
	with uigen.group(ui, 'Advanced'):
		ui.parameter('glintRolloff', 'float', 'Glint rolloff', connectible=False)
		ui.parameter('transmissionRolloff', 'float', 'Transmission rolloff', connectible=False)
		ui.parameter('singleSaturation', 'float', 'Highlight saturation')
		ui.parameter('multipleSaturation', 'float', 'Diffuse saturation')

with uigen.group(ui, 'Diffuse', False):
	ui.parameter('diffuseStrength', 'float', 'Strength')
	ui.parameter('diffuseColor', 'rgb', 'Color')
	ui.parameter('diffuseForward', 'float', 'Forward scattering')
	ui.parameter('diffuseBack', 'float', 'Back scattering')
	with uigen.group(ui, 'Advanced'):
		ui.parameter('dualDepth', 'int', 'Dual depth')

with uigen.group(ui, 'Specular 1', False):
	ui.parameter('specular1Strength', 'float', 'Strength')
	ui.parameter('specular1Color', 'rgb', 'Tint')
	ui.parameter('specular1WidthScale', 'float', 'Width scale')

with uigen.group(ui, 'Specular 2', False):
	ui.parameter('specular2Strength', 'float', 'Strength')
	ui.parameter('specular2Color', 'rgb', 'Tint')
	ui.parameter('specular2WidthScale', 'float', 'Width scale')
	with uigen.group(ui, 'Glints', False):
		ui.parameter('glintStrength', 'float', 'Strength')
		ui.parameter('glintTexture', 'float', 'Texture')
		ui.parameter('twist', 'float', 'Twist')

with uigen.group(ui, 'Transmission', False):
	ui.parameter('transmissionStrength', 'float', 'Strength')
	ui.parameter('transmissionColor', 'rgb', 'Tint')
	ui.parameter('transmissionWidthScale', 'float', 'Width scale')

with uigen.group(ui, 'AOVs'):
	ui.aov('aov_diffuse_color', 'rgb', 'Diffuse color')
	ui.aov('aov_direct_diffuse', 'rgb', 'Direct diffuse')
	ui.aov('aov_indirect_diffuse', 'rgb', 'Indirect diffuse')
	ui.aov('aov_direct_local', 'rgb', 'Direct local')
	ui.aov('aov_indirect_local', 'rgb', 'Indirect local')
	ui.aov('aov_direct_global', 'rgb', 'Direct global')
	ui.aov('aov_indirect_global', 'rgb', 'Indirect global')
	ui.aov('aov_direct_specular', 'rgb', 'Direct specular')
	ui.aov('aov_indirect_specular', 'rgb', 'Indirect specular')
	ui.aov('aov_direct_specular_2', 'rgb', 'Direct specular 2')
	ui.aov('aov_indirect_specular_2', 'rgb', 'Indirect specular 2')
	ui.aov('aov_direct_transmission', 'rgb', 'Direct transmission')
	ui.aov('aov_indirect_transmission', 'rgb', 'Indirect transmission')
	ui.aov('aov_light_group_1', 'rgb', 'Light group [1]')
	ui.aov('aov_light_group_2', 'rgb', 'Light group [2]')
	ui.aov('aov_light_group_3', 'rgb', 'Light group [3]')
	ui.aov('aov_light_group_4', 'rgb', 'Light group [4]')
	ui.aov('aov_light_group_5', 'rgb', 'Light group [5]')
	ui.aov('aov_light_group_6', 'rgb', 'Light group [6]')
	ui.aov('aov_light_group_7', 'rgb', 'Light group [7]')
	ui.aov('aov_light_group_8', 'rgb', 'Light group [8]')
	ui.aov('aov_id_1', 'rgb', 'ID [1]')
	ui.aov('aov_id_2', 'rgb', 'ID [2]')
	ui.aov('aov_id_3', 'rgb', 'ID [3]')
	ui.aov('aov_id_4', 'rgb', 'ID [4]')
	ui.aov('aov_id_5', 'rgb', 'ID [5]')
	ui.aov('aov_id_6', 'rgb', 'ID [6]')
	ui.aov('aov_id_7', 'rgb', 'ID [7]')
	ui.aov('aov_id_8', 'rgb', 'ID [8]')

