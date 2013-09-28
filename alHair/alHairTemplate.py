import pymel.core as pm
from alShaders import alShadersTemplate

class AEalHairTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.beginLayout("Matte", collapse=True)
        self.addControl("aiEnableMatte", label="Enable Matte")
        self.addControl("aiMatteColor", label="Matte Color")
        self.addControl("aiMatteColorA", label="Matte Opacity")
        self.endLayout()

        self.beginLayout('Fibre properties', collapse=False)
        self.addControl('hairColor', label='Color')
        self.addControl('specularWidth', label='Highlight width')
        self.addControl('specularShift', label='Highlight shift')
        self.addControl('opacity')
        self.beginLayout('Advanced')
        self.addControl('transmissionRolloff', label='Transmission rolloff')
        self.addControl('glintRolloff', label='Glint rolloff')
        self.addControl('singleSaturation')
        self.addControl('multipleSaturation')
        self.endLayout()
        self.endLayout() # end Fibre properties

        self.beginLayout('Diffuse', collapse=False)
        self.addControl('diffuseStrength', label='Strength')
        self.addControl('diffuseColor', label='Tint')
        self.addControl('diffuseForward', label='Forward scattering')
        self.addControl('diffuseBack', label="Back scattering")

        self.beginLayout('Advanced')
        self.addControl('dualDepth')
        self.endLayout()
        self.endLayout() # end Dual scattering

        self.beginLayout('Specular 1', collapse=False)
        self.addControl('specular1Strength', label='Strength')
        self.addControl('specular1Color', label='Tint')
        self.addControl('specular1WidthScale', label='Width scale')
        self.endLayout() # end Specular 1

        self.beginLayout('Specular 2', collapse=False)
        self.addControl('specular2Strength', label='Strength')
        self.addControl('specular2Color', label='Tint')
        self.addControl('specular2WidthScale', label='Width scale')
        self.beginLayout('Glints')
        self.addControl('glintStrength', label='Strength')
        self.addControl('glintTexture', label='Texture')
        self.addControl('twist', label='Twist')
        self.endLayout()
        self.beginLayout('Advanced')
        self.endLayout()
        self.endLayout() # end Specular 2

        self.beginLayout('Transmission', collapse=False)
        self.addControl('transmissionStrength', label='Strength')
        self.addControl('transmissionColor', label='Tint')
        self.addControl('transmissionWidthScale', label='Width scale')
        self.beginLayout('Advanced')
        
        self.endLayout()
        self.endLayout() # end Transmission

        self.beginLayout('AOVs', collapse=True)
        #self.addControl('lightGroupsIndirect', label='Indirect light groups')
        self.addControl('aov_diffuse_color', label='Diffuse color')
        self.addControl('aov_direct_diffuse', label='Direct diffuse')
        self.addControl('aov_indirect_diffuse', label='Indirect diffuse')
        self.addControl('aov_direct_local', label='Direct local')
        self.addControl('aov_indirect_local', label='Indirect local')
        self.addControl('aov_direct_global', label='Direct global')
        self.addControl('aov_indirect_global', label='Indirect global')
        self.addControl('aov_direct_specular', label='Direct specular')
        self.addControl('aov_indirect_specular', label='Indirect specular')
        self.addControl('aov_direct_specular_2', label='Direct specular 2')
        self.addControl('aov_indirect_specular_2', label='Indirect specular 2')
        self.addControl('aov_direct_glint', label='Direct glint')
        self.addControl('aov_indirect_glint', label='Indirect glint')
        self.addControl('aov_direct_transmission', label='Direct transmission')
        self.addControl('aov_indirect_transmission', label='Indirect transmission')
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

        self.beginLayout('Advanced', collapse=True)
        self.addControl('extraSamples', label='Extra samples')
        self.addControl('MIS')
        self.endLayout() # end Advacned

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()