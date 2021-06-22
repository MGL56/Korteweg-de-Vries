# Korteweg-de-Vries
Résolution de l'équation Korteweg-de-Vries

On considère l'équation KdV sur un domaine [-L, L]:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\zeta}{\partial&space;t}&space;&plus;&space;\frac{\partial&space;\zeta}{\partial&space;x}&space;&plus;\frac{3\epsilon}{2}\frac{\partial&space;\zeta}{\partial&space;x}\zeta&space;&plus;&space;\frac{\epsilon}{6}\frac{\partial^3&space;\zeta}{\partial&space;x^3}&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\zeta}{\partial&space;t}&space;&plus;&space;\frac{\partial&space;\zeta}{\partial&space;x}&space;&plus;\frac{3\epsilon}{2}\frac{\partial&space;\zeta}{\partial&space;x}\zeta&space;&plus;&space;\frac{\epsilon}{6}\frac{\partial^3&space;\zeta}{\partial&space;x^3}&space;=&space;0" title="\frac{\partial \zeta}{\partial t} + \frac{\partial \zeta}{\partial x} +\frac{3\epsilon}{2}\frac{\partial \zeta}{\partial x}\zeta + \frac{\epsilon}{6}\frac{\partial^3 \zeta}{\partial x^3} = 0" /></a>

On discrétise cette équation avec des différences finies, l'inconnue est approximée sur les points

<a href="https://www.codecogs.com/eqnedit.php?latex=x{_{i}}&space;=&space;-L&space;&plus;&space;2\frac{i}{N}L,&space;i&space;=&space;0..N-1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x{_{i}}&space;=&space;-L&space;&plus;&space;2\frac{i}{N}L,&space;i&space;=&space;0..N-1" title="x{_{i}} = -L + 2\frac{i}{N}L, i = 0..N-1" /></a>
