import pymel.core as pm
from mtoa.ui.ae.shaderTemplate import ShaderAETemplate

UI = {
    'float': pm.attrFieldSliderGrp,
    'rgb': pm.attrColorSliderGrp,
}

class Param:
    name = ''
    label = ''
    annotation = ''
    ptype = ''
    presets = None
    precision = 4

    def __init__(self, n, l, a, p, presets=None, precision=4):
        self.name = n
        self.label = l
        self.annotation = a
        self.ptype = p
        self.presets = presets
        self.precision = precision

def setPresetFlt(ctrl, value):
    attr = pm.attrFieldSliderGrp(ctrl, query=True, attribute=True)
    pm.setAttr(attr, value)

def setPresetRgb(ctrl, value):
    attr = pm.attrColorSliderGrp(ctrl, query=True, attribute=True)
    pm.setAttr(attr + 'R', value[0])
    pm.setAttr(attr + 'G', value[1])
    pm.setAttr(attr + 'B', value[2])

class alShadersTemplate(ShaderAETemplate):

    def customCreateFlt(self, attr):
        pname = attr.split('.')[-1]
        ptype = self.params[pname].ptype
        plabel = self.params[pname].label
        pann = self.params[pname].annotation
        presets = self.params[pname].presets
        precision = self.params[pname].precision
        self.controls[pname] = pm.attrFieldSliderGrp(pname + 'Ctrl', attribute=attr, label=plabel, annotation=pann, precision=precision)
        if presets is not None:
            pm.attrFieldSliderGrp(self.controls[pname], edit=True, bgc=(.22, .22, .22))
            pm.popupMenu()
            for k in sorted(presets, key=presets.get):
                pm.menuItem(label=k, command=pm.Callback(setPresetRgb, self.controls[pname], presets[k]))


    def customUpdateFlt(self, attr):
        pname = attr.split('.')[-1]
        ptype = self.params[pname].ptype
        pm.attrFieldSliderGrp(self.controls[pname], edit=True, attribute=attr)

    def customCreateRgb(self, attr):
        pname = attr.split('.')[-1]
        ptype = self.params[pname].ptype
        plabel = self.params[pname].label
        pann = self.params[pname].annotation
        presets = self.params[pname].presets
        self.controls[pname] = pm.attrColorSliderGrp(pname + 'Ctrl', attribute=attr, label=plabel, annotation=pann)
        if presets is not None:
            pm.attrColorSliderGrp(self.controls[pname], edit=True, bgc=(.22, .22, .22))
            pm.popupMenu()
            for k in sorted(presets, key=presets.get):
                pm.menuItem(label=k, command=pm.Callback(setPresetRgb, self.controls[pname], presets[k]))

    def customUpdateRgb(self, attr):
        pname = attr.split('.')[-1]
        ptype = self.params[pname].ptype
        pm.attrColorSliderGrp(self.controls[pname], edit=True, attribute=attr)


    def addCustomFlt(self, param):
        self.addCustom(param, self.customCreateFlt, self.customUpdateFlt)

    def addCustomRgb(self, param):
        self.addCustom(param, self.customCreateRgb, self.customUpdateRgb)

    def addRemapControls(self):
        self.beginLayout('Remap', collapse=True)
        self.addControl('RMPinputMin', label='Input min')
        self.addControl('RMPinputMax', label='Input max')

        self.beginLayout('Contrast', collapse=False)
        self.addControl('RMPcontrast', label='Contrast')
        self.addControl('RMPcontrastPivot', label='Pivot')
        self.endLayout() # end Contrast

        self.beginLayout('Bias and gain', collapse=False)
        self.addControl('RMPbias', label='Bias')
        self.addControl('RMPgain', label='Gain')
        self.endLayout() # end Bias and gain

        self.addControl('RMPoutputMin', label='Output min')
        self.addControl('RMPoutputMax', label='Output max')

        self.beginLayout('Clamp', collapse=False)
        self.addControl('RMPclampEnable', label='Enable')
        self.addControl('RMPthreshold', label='Expand')
        self.addControl('RMPclampMin', label='Min')
        self.addControl('RMPclampMax', label='Max')
        self.endLayout() # end Clamp

        self.endLayout() # end Remap

    
    