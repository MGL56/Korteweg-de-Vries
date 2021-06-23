# Korteweg-de-Vries
Résolution de l'équation Korteweg-de-Vries

On considère l'équation KdV sur un domaine [-L, L]:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\zeta}{\partial&space;t}&space;&plus;&space;\frac{\partial&space;\zeta}{\partial&space;x}&space;&plus;\frac{3\epsilon}{2}\frac{\partial&space;\zeta}{\partial&space;x}\zeta&space;&plus;&space;\frac{\epsilon}{6}\frac{\partial^3&space;\zeta}{\partial&space;x^3}&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\zeta}{\partial&space;t}&space;&plus;&space;\frac{\partial&space;\zeta}{\partial&space;x}&space;&plus;\frac{3\epsilon}{2}\frac{\partial&space;\zeta}{\partial&space;x}\zeta&space;&plus;&space;\frac{\epsilon}{6}\frac{\partial^3&space;\zeta}{\partial&space;x^3}&space;=&space;0" title="\frac{\partial \zeta}{\partial t} + \frac{\partial \zeta}{\partial x} +\frac{3\epsilon}{2}\frac{\partial \zeta}{\partial x}\zeta + \frac{\epsilon}{6}\frac{\partial^3 \zeta}{\partial x^3} = 0" /></a>

On discrétise cette équation avec des différences finies, l'inconnue est approximée sur les points

<a href="https://www.codecogs.com/eqnedit.php?latex=x{_{i}}&space;=&space;-L&space;&plus;&space;2\frac{i}{N}L,&space;i&space;=&space;0..N-1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x{_{i}}&space;=&space;-L&space;&plus;&space;2\frac{i}{N}L,&space;i&space;=&space;0..N-1" title="x{_{i}} = -L + 2\frac{i}{N}L, i = 0..N-1" /></a>

On remarquera que l'extémité L n'est pas discrétisée puisque on sait que la solution à ce point est égale à la solution au point -L. On introduit les deux matrices D1 et D3 qui discrétisent la dérivée première et troisième avec un simple schéma centré d'ordre 2:

<a href="https://www.codecogs.com/eqnedit.php?latex=D_{1}U&space;=&space;\frac{U_{i&plus;1}-U_{i-1}}{2\Delta&space;x}&space;\\" target="_blank"><img src="https://latex.codecogs.com/gif.latex?D_{1}U&space;=&space;\frac{U_{i&plus;1}-U_{i-1}}{2\Delta&space;x}&space;\\" title="D_{1}U = \frac{U_{i+1}-U_{i-1}}{2\Delta x} \\" /></a>


<a href="https://www.codecogs.com/eqnedit.php?latex=D_{3}U&space;=&space;\frac{U_{i&plus;2}-2U_{i&plus;1}&space;&plus;&space;2U_{i-1}&space;-U_{i-2}}{2\Delta&space;x^{3}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?D_{3}U&space;=&space;\frac{U_{i&plus;2}-2U_{i&plus;1}&space;&plus;&space;2U_{i-1}&space;-U_{i-2}}{2\Delta&space;x^{3}}" title="D_{3}U = \frac{U_{i+2}-2U_{i+1} + 2U_{i-1} -U_{i-2}}{2\Delta x^{3}}" /></a>

Du fait des conditions périodiques, les matrices D1 et D3 sont égales à:

<a href="https://www.codecogs.com/eqnedit.php?latex=D1&space;=&space;\frac{1}{2\Delta&space;x}&space;\begin{pmatrix}&space;0&1&0&\cdots&0&-1\\&space;-1&&space;0&&space;1&&space;\ddots&space;&&space;0&space;&0&space;\\&space;0&&space;-1&&space;0&&space;\ddots&space;&&space;\vdots&space;&0&space;\\&space;\vdots&space;&\ddots&space;&&space;\ddots&space;&\ddots&space;&\ddots&space;&\vdots&space;\\&space;0&space;&&space;0&&space;\cdots&space;&&space;-1&0&space;&&space;1\\&space;1&space;&0&space;&\cdots&space;&&space;0&space;&-1&space;&&space;0&space;\end{pmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?D1&space;=&space;\frac{1}{2\Delta&space;x}&space;\begin{pmatrix}&space;0&1&0&\cdots&0&-1\\&space;-1&&space;0&&space;1&&space;\ddots&space;&&space;0&space;&0&space;\\&space;0&&space;-1&&space;0&&space;\ddots&space;&&space;\vdots&space;&0&space;\\&space;\vdots&space;&\ddots&space;&&space;\ddots&space;&\ddots&space;&\ddots&space;&\vdots&space;\\&space;0&space;&&space;0&&space;\cdots&space;&&space;-1&0&space;&&space;1\\&space;1&space;&0&space;&\cdots&space;&&space;0&space;&-1&space;&&space;0&space;\end{pmatrix}" title="D1 = \frac{1}{2\Delta x} \begin{pmatrix} 0&1&0&\cdots&0&-1\\ -1& 0& 1& \ddots & 0 &0 \\ 0& -1& 0& \ddots & \vdots &0 \\ \vdots &\ddots & \ddots &\ddots &\ddots &\vdots \\ 0 & 0& \cdots & -1&0 & 1\\ 1 &0 &\cdots & 0 &-1 & 0 \end{pmatrix}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=D3=\frac{1}{2\Delta&space;x^{3}}&space;\begin{pmatrix}&space;0&-2&1&0&\cdots&-1&2\\&space;2&&space;0&&space;-2&&space;1&&space;\ddots&space;&&space;0&space;&-&space;\\&space;-1&&space;2&&space;0&&space;\ddots&space;&&space;\ddots&space;&&space;\vdots&space;&0&space;\\&space;0&\ddots&space;&&space;\ddots&space;&\ddots&space;&\ddots&space;&\vdots&space;&\vdots\\&space;0&space;&&space;\ddots&space;&&space;-1&2&0&-2&space;&&space;-1\\&space;1&space;&0&space;&\cdots&space;&&space;-1&space;&2&space;&&space;0&space;&-2&space;\\&space;-2&1&0&\cdots&-1&2&0&space;\end{pmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?D3=\frac{1}{2\Delta&space;x^{3}}&space;\begin{pmatrix}&space;0&-2&1&0&\cdots&-1&2\\&space;2&&space;0&&space;-2&&space;1&&space;\ddots&space;&&space;0&space;&-&space;\\&space;-1&&space;2&&space;0&&space;\ddots&space;&&space;\ddots&space;&&space;\vdots&space;&0&space;\\&space;0&\ddots&space;&&space;\ddots&space;&\ddots&space;&\ddots&space;&\vdots&space;&\vdots\\&space;0&space;&&space;\ddots&space;&&space;-1&2&0&-2&space;&&space;-1\\&space;1&space;&0&space;&\cdots&space;&&space;-1&space;&2&space;&&space;0&space;&-2&space;\\&space;-2&1&0&\cdots&-1&2&0&space;\end{pmatrix}" title="D3=\frac{1}{2\Delta x^{3}} \begin{pmatrix} 0&-2&1&0&\cdots&-1&2\\ 2& 0& -2& 1& \ddots & 0 &- \\ -1& 2& 0& \ddots & \ddots & \vdots &0 \\ 0&\ddots & \ddots &\ddots &\ddots &\vdots &\vdots\\ 0 & \ddots & -1&2&0&-2 & -1\\ 1 &0 &\cdots & -1 &2 & 0 &-2 \\ -2&1&0&\cdots&-1&2&0 \end{pmatrix}" /></a>
