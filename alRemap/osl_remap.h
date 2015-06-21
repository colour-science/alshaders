#define REMAP_FLOAT_DECLARE_PARAMS \
float RMPinputMin = 0         [[string label = "Input min", string page="Remap"]], \
float RMPinputMax = 1         [[string label = "Input max", string page="Remap"]], \
float RMPcontrast = 1      [[string label = "Contrast", string page="Remap/Contrast"]], \
float RMPcontrastPivot = 0.5  [[string label = "Pivot", string page="Remap/Contrast"]], \
float RMPbias = 0.5           [[string label = "Bias", string page="Remap/Bias and Gain"]], \
float RMPgain = 0.5           [[string label = "Gain", string page="Remap/Bias and Gain"]], \
float RMPoutputMin = 0        [[string label = "Output min", string page="Remap"]], \
float RMPoutputMax = 1        [[string label = "Output max", string page="Remap"]], \
int RMPclampEnable = 1        [[string label = "Clamp", string widget = "checkBox", string page="Remap/Clamp"]], \
int RMPthreshold = 0             [[string label = "Expand", string widget = "checkBox", string page="Remap/Clamp"]], \
float RMPclampMin = 0         [[string label = "Min", string page="Remap/Clamp"]], \
float RMPclampMax = 1         [[string label = "Max", string page="Remap/Clamp"]]

struct RemapFloatParams
{
   float inputMin;
   float inputMax;
   float contrastVal;
   float contrastPivot;
   float bias;
   float gain;
   float outputMin;
   float outputMax;
   int clampEnable;
   int expand;
   float clampMin;
   float clampMax;
};

#define REMAP_FLOAT_CREATE \
{RMPinputMin, RMPinputMax, RMPcontrast, RMPcontrastPivot, RMPbias, RMPgain, RMPoutputMin, RMPoutputMax, RMPclampEnable, RMPthreshold, RMPclampMin, RMPclampMax}

float contrast(float input, float c, float pivot)
{
    if (c == 1.0) return input;

    return (input-pivot)*c + pivot;
}

float bias(float f, float b)
{
    if (b > 0.0) return pow(f, log(b)/log(0.5));
    else return 0.0;
}

float biasandgain(float f, float b, float g)
{
    if (f < 0) return f;

    if (b != 0.5)
    {
        f = bias(f, b);
    }
    if (g != 0.5)
    {
        if (f < 0.5) f = 0.5 * bias(2.0*f, 1.0-g);
        else f = 1.0 - bias(2.0 - 2.0*f, 1.0-g)*0.5;
    }
    return f;
}

float lerp(float a, float b, float t)
{
   return (1-t)*a + t*b;
}

float remap(float input, RemapFloatParams p)
{
   float f = (input-p.inputMin)/(p.inputMax-p.inputMin);
   f = contrast(f, p.contrastVal, p.contrastPivot);
   f = biasandgain(f, p.bias, p.gain);
   f = lerp(p.outputMin, p.outputMax, f);
   if (p.clampEnable)
   {
      f = min(p.clampMax, f);
      f = max(p.clampMin, f);
      if (p.expand)
      {
         f = (f-p.clampMin)/(p.clampMax-p.clampMin);
      }
   }
   return f;
}
