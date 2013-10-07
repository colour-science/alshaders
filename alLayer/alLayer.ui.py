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
	'soft_version':1
})

ui.parameter('layer1', 'rgb', 'Layer 1')
ui.parameter('layer2', 'rgb', 'Layer 2')
ui.parameter('mix', 'float', 'Mix')
ui.parameter('debug', 'enum', 'Debug', enum_names=[
	"off",
	"layer1",
	"layer2",
	"mixer"
])

with uigen.group(ui, 'AOVs'):
	ui.parameter('standardCompatibleAOVs', 'bool', 'Write standard AOVs only')
	ui.parameter('aov_diffuse_color', 'rgb', 'Diffuse color')
	ui.parameter('aov_direct_diffuse', 'rgb', 'Direct diffuse')
	ui.parameter('aov_direct_diffuse_raw', 'rgb', 'Direct diffuse (raw)')
	ui.parameter('aov_indirect_diffuse', 'rgb', 'Indirect diffuse')
	ui.parameter('aov_indirect_diffuse_raw', 'rgb', 'Indirect diffuse (raw)')
	ui.parameter('aov_direct_backlight', 'rgb', 'Direct backlight')
	ui.parameter('aov_indirect_backlight', 'rgb', 'Indirect backlight')
	ui.parameter('aov_direct_specular', 'rgb', 'Direct specular')
	ui.parameter('aov_indirect_specular', 'rgb', 'Indirect specular')
	ui.parameter('aov_direct_specular_2', 'rgb', 'Direct specular 2')
	ui.parameter('aov_indirect_specular_2', 'rgb', 'Indirect specular 2')
	ui.parameter('aov_single_scatter', 'rgb', 'Single scatter')
	ui.parameter('aov_sss', 'rgb', 'SSS')
	ui.parameter('aov_refraction', 'rgb', 'Refraction')
	ui.parameter('aov_emission', 'rgb', 'Emission')
	ui.parameter('aov_uv', 'rgb', 'UV')
	ui.parameter('aov_depth', 'rgb', 'Depth')
	ui.parameter('aov_light_group_1', 'rgb', 'Light group [1]')
	ui.parameter('aov_light_group_2', 'rgb', 'Light group [2]')
	ui.parameter('aov_light_group_3', 'rgb', 'Light group [3]')
	ui.parameter('aov_light_group_4', 'rgb', 'Light group [4]')
	ui.parameter('aov_light_group_5', 'rgb', 'Light group [5]')
	ui.parameter('aov_light_group_6', 'rgb', 'Light group [6]')
	ui.parameter('aov_light_group_7', 'rgb', 'Light group [7]')
	ui.parameter('aov_light_group_8', 'rgb', 'Light group [8]')
	ui.parameter('aov_id_1', 'rgb', 'ID [1]')
	ui.parameter('aov_id_2', 'rgb', 'ID [2]')
	ui.parameter('aov_id_3', 'rgb', 'ID [3]')
	ui.parameter('aov_id_4', 'rgb', 'ID [4]')
	ui.parameter('aov_id_5', 'rgb', 'ID [5]')
	ui.parameter('aov_id_6', 'rgb', 'ID [6]')
	ui.parameter('aov_id_7', 'rgb', 'ID [7]')
	ui.parameter('aov_id_8', 'rgb', 'ID [8]')
