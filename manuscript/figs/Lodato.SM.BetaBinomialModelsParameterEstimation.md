# Model parameters from Lodato Supplementary

extracted from positional values in Inkscape

## Fig S5C, beta-binomial mixture model with linear scaling

### model

$P_{het}(alt.read.depth = m|tot.depth = n) = w_n*BB(m;\alpha = \alpha_1(n); \beta = \alpha_1(n), n = n) + (1 - w_n)*BB(m; \alpha = \alpha_2(n), \beta = \alpha_2(n), n = n)$

Then, we assume a linear trend of the mixture model parameters with respect to n.
$\alpha_1(n) = a_1*n + b_1$
$\alpha_2(n) = a_2*n + b_2$
$w(n) = a_3*n + b_3$

### parameter value scale:

(124,210-33,692)/4 = 22,6295 px per unit

60 - 5 = 55 units on total depth scale

### alpha1

slope: ((72,145-(64,77-64,04))/22,6295)/55 = 0,057378844 units per total depth

(heightOfLine-(TopOfLineAtDepth10-BottomOfLineAtDepth10))/ParameterScaleValue

intercept: (55,340-33,692)/22,6295 - 0,057378844*5 = 0,669733191 units

(YAtDepth5-YOf0ParameterValueTick)/ParameterScaleValue - Slope*Depth5

### alpha2

slope: ((4,635-(44,2-43,59))/22,6295)/55 = 0,003233912 units per total depth

intercept: (43,093-33,692)/22,6295 - 0,003233912*5 = 0,399261625 units

### w

slope: ((1,303-(46,68-46,06))/22,6295)/55 = 0,000548761 units per total depth

intercept: (45,983-33,692)/22,6295 - 0,000548761*5 = 0,540396786 units


## Fig S5D, beta-binomial model

### model

$P_{nonvariant}(alt.read.depth = m|tot.depth = n) = BB(m;\alpha(n),\beta(n),n)

$\alpha(n) = a_1*n + b_1$
$\beta(n) = a_2*n + b_2$

### parameter value scale:

(127,583-33,942)/4 = 23,41025 px per unit

### alpha

slope: ((-0,505+(36,01-35,54))/23,41025)/55 = -0,000027183 units per total depth

intercept: (35,509-33,942)/23,41025 - (-0,000027183)*60 = 0,068567471 units

### beta

slope: ((10,228-(92,04-91,41))/23,41025)/55 = 0,007454388 units per total depth

intercept: (90,238-33,942)/23,41025 - 0,007454388*5 = 2,367486659 units
