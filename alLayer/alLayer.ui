import uigen

ui.shader({
	'name':'alLayer',
	'description':'Layer two shaders together',
	'output':'rgb',
	'maya_name':'alLayer',
	'maya_classification':'shader/surface',
	'maya_id':'0x00116407',
	'maya_swatch':True,
	'maya_matte':True,
	'maya_bump':False,
	'soft_name':'ALS_Layer',
	'soft_classification':'material',
	'soft_version':1,
	'help_url':'https://bitbucket.org/anderslanglands/alshaders/wiki/alLayer'
})

ui.parameter('layer1', 'rgb', (0.0, 0.0, 0.0), 'Layer 1')
ui.parameter('layer2', 'rgb', (0.0, 0.0, 0.0), 'Layer 2')
ui.parameter('mix', 'float', 0.0, 'Mix')
ui.parameter('debug', 'enum', 'off', 'Debug', enum_names=[
	"off",
	"layer1",
	"layer2",
	"mixer"
])

with uigen.group(ui, 'AOVs'):
	ui.parameter('standardCompatibleAOVs', 'bool', False, 'Write standard AOVs only')
	ui.aov('aov_diffuse_color', 'rgb', 'Diffuse color')
	ui.aov('aov_direct_diffuse', 'rgb', 'Direct diffuse')
	ui.aov('aov_direct_diffuse_raw', 'rgb', 'Direct diffuse (raw)')
	ui.aov('aov_indirect_diffuse', 'rgb', 'Indirect diffuse')
	ui.aov('aov_indirect_diffuse_raw', 'rgb', 'Indirect diffuse (raw)')
	ui.aov('aov_direct_backlight', 'rgb', 'Direct backlight')
	ui.aov('aov_indirect_backlight', 'rgb', 'Indirect backlight')
	ui.aov('aov_direct_specular', 'rgb', 'Direct specular')
	ui.aov('aov_indirect_specular', 'rgb', 'Indirect specular')
	ui.aov('aov_direct_specular_2', 'rgb', 'Direct specular 2')
	ui.aov('aov_indirect_specular_2', 'rgb', 'Indirect specular 2')
	ui.aov('aov_single_scatter', 'rgb', 'Single scatter')
	ui.aov('aov_sss', 'rgb', 'SSS')
	ui.aov('aov_refraction', 'rgb', 'Refraction')
	ui.aov('aov_emission', 'rgb', 'Emission')
	ui.aov('aov_uv', 'rgb', 'UV')
	ui.aov('aov_depth', 'rgb', 'Depth')
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
