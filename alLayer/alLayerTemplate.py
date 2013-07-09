import pymel.core as pm
from alShaders import alShadersTemplate

class AEalLayerTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('layer1', label='Layer 1')
        self.addControl('layer2', label='Layer 2')
        self.addControl('mix', label='Mix')
        self.addControl('debug', label='Debug')

        self.beginLayout('AOVs', collapse=True)
        self.addControl('standardCompatibleAOVs', label='Write standard AOVs only')
        self.addControl('aov_diffuse_color', label='Diffuse color')
        self.addControl('aov_direct_diffuse', label='Direct diffuse')
        self.addControl('aov_direct_diffuse_raw', label='Direct diffuse (raw)')
        self.addControl('aov_indirect_diffuse', label='Indirect diffuse')
        self.addControl('aov_indirect_diffuse_raw', label='Indirect diffuse (raw)')
        self.addControl('aov_direct_backlight', label='Direct backlight')
        self.addControl('aov_indirect_backlight', label='Indirect backlight')
        self.addControl('aov_direct_specular', label='Direct specular')
        self.addControl('aov_indirect_specular', label='Indirect specular')
        self.addControl('aov_direct_specular_2', label='Direct specular 2')
        self.addControl('aov_indirect_specular_2', label='Indirect specular 2')
        self.addControl('aov_single_scatter', label='Single scatter')
        self.addControl('aov_sss', label='SSS')
        self.addControl('aov_refraction', label='Refraction')
        self.addControl('aov_emission', label='Emission')
        self.addControl('aov_uv', label='UV')
        self.addControl('aov_depth', label='Depth')
        self.addControl('aov_light_group_1', label='Light group [1]')
        self.addControl('aov_light_group_2', label='Light group [2]')
        self.addControl('aov_light_group_3', label='Light group [3]')
        self.addControl('aov_light_group_4', label='Light group [4]')
        self.addControl('aov_light_group_5', label='Light group [5]')
        self.addControl('aov_light_group_6', label='Light group [6]')
        self.addControl('aov_light_group_7', label='Light group [7]')
        self.addControl('aov_light_group_8', label='Light group [8]')
        self.addControl('aov_id_1', label='ID [1]')
        self.addControl('aov_id_2', label='ID [2]')
        self.addControl('aov_id_3', label='ID [3]')
        self.addControl('aov_id_4', label='ID [4]')
        self.addControl('aov_id_5', label='ID [5]')
        self.addControl('aov_id_6', label='ID [6]')
        self.addControl('aov_id_7', label='ID [7]')
        self.addControl('aov_id_8', label='ID [8]')
        self.endLayout() # end AOVs

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()