<map version="0.9.0">
<!--To view this file, download free mind mapping software Freeplane from http://freeplane.sourceforge.net -->
<node TEXT="lKTA" ID="ID_1723255651" CREATED="1283093380553" MODIFIED="1456600614781">
<hook NAME="MapStyle" max_node_width="600"/>
<hook NAME="AutomaticEdgeColor" COUNTER="4"/>
<node TEXT="Analisis de resolucion geneX&#xa; + estructura de particion&#xa; + BHI (solo BHI vs Size)&#xa;Comparar particiones usando ecdf y compare partitions?" POSITION="right" ID="ID_1561909303" CREATED="1456600301595" MODIFIED="1456684636904">
<edge COLOR="#00ff00"/>
<icon BUILTIN="help"/>
<node TEXT="" ID="ID_715396996" CREATED="1456600806947" MODIFIED="1456600903566">
<node TEXT="kmeans (nro variable de clusters, nro sugerido por metricas)" ID="ID_784775548" CREATED="1456600845188" MODIFIED="1456600866854">
<node TEXT="distancia con coordenadas estandarizadas" ID="ID_540624027" CREATED="1456777227264" MODIFIED="1456777242236"/>
<node TEXT="barrido en k , figura de merito, k=2" ID="ID_1400447344" CREATED="1456777032146" MODIFIED="1456777050482"/>
<node TEXT="problema: estructuras muy gdes (rho?), dificil interpretacion" ID="ID_1585732741" CREATED="1456777051392" MODIFIED="1456777265276"/>
</node>
<node TEXT="cutreedynamic" ID="ID_1425066013" CREATED="1456600868057" MODIFIED="1456600879236">
<node TEXT="DS4" ID="ID_734505956" CREATED="1456600967074" MODIFIED="1456600969232"/>
<node TEXT="DS1" ID="ID_31693319" CREATED="1456600969804" MODIFIED="1456600975854"/>
</node>
<node TEXT="plot ecdf y compare_partitions comparando los 3&#xa;ver que mas o menos DS4 C DS1 C kmeans" ID="ID_1516495765" CREATED="1456777473128" MODIFIED="1456781519937" STYLE="bubble"/>
<node TEXT="conclusion: con el espacio geneX podemos encontrar estructura que no siempre es facil de interpretar ....por un lado esto sucede cuando hay bodoques, pero tampoco es claro encontrar la resolucion optima de analisis.&#xa;en el capitulo siguiente introduciremos una herramienta que busca cuantificar la homogeneidad bioplogica de particiones para encontrar dicha escala." ID="ID_1227325308" CREATED="1456600983667" MODIFIED="1456777893177"/>
<node TEXT="" ID="ID_1588634234" CREATED="1456777776072" MODIFIED="1456777776072"/>
</node>
</node>
<node TEXT="analisis de BHI como medida (usando algun/os tratamientos DS1)" POSITION="right" ID="ID_1128389210" CREATED="1456599998239" MODIFIED="1456600783855" VSHIFT="-100">
<edge COLOR="#ff0000"/>
<font BOLD="false"/>
<node TEXT="variantes de BHI:&#xa;BHI-delta, BHI_IC, BHI_Resnick&#xa;(estos ser&#xed;an los que llamamos BHI, BHIpesado, BHIresnick?)" ID="ID_965197812" CREATED="1456600658698" MODIFIED="1456683043235">
<icon BUILTIN="help"/>
</node>
<node TEXT="analisis de correlacion...chequear si se justifica usar BHI simple&#xa;(esto es ver que todos los tipos de BHI tienen la misma dependencia con el tama&#xf1;o de los clusters y que por eso vamos a usar el m&#xe1;s sencillo?)" ID="ID_657151194" CREATED="1456600681708" MODIFIED="1456781544089" STYLE="bubble">
<icon BUILTIN="help"/>
</node>
<node TEXT="hablar de controles nulos" ID="ID_1915119892" CREATED="1456600717764" MODIFIED="1456600727901">
<node TEXT="independent clusters" ID="ID_128615720" CREATED="1456600727903" MODIFIED="1456600734445">
<node TEXT="BHINullGenerator.R&#xa;tomando 6000 genes que se movieron alguna vez&#xa;1000-bootstrap: &lt;BHI&gt;~cte, sdBHI=sdBHI(size)&#xa;fiteado por 2 powerlaws:[1,50][51,...]." ID="ID_816359629" CREATED="1456605012755" MODIFIED="1456605181403"/>
</node>
<node TEXT="label random shuffling" ID="ID_965189990" CREATED="1456600735223" MODIFIED="1456600757983">
<node TEXT="definir. Comparar con la otra NULL," ID="ID_264148802" CREATED="1456605184017" MODIFIED="1456605260928"/>
</node>
</node>
<node TEXT="analisis de bhi para kmeans, ds1, ds4&#xa;tres paneles  bhi vs size para los clusters detectados por kmeans, ds1 y ds4 respectivamente&#xa;con colores por tratamiento" ID="ID_18375450" CREATED="1456777894944" MODIFIED="1456781559721" STYLE="bubble"/>
<node TEXT="deberiamos poder decir que en general no es buno el valor de bhi&#xa;y que mayor resolucion en geneX no mejora.&#xa;cuantificar con un observable para ver esto&#xa;&#xa;conclusion: con el espacio geneX podemos encontrar estructura que no siempre es facil de interpretar a la luz de conocimiento almacenado en estructura de GO....por un lado esto sucede cuando hay bodoques, pero tampoco mejora demasiado cuando se consideran particiones con estructuras de mas resolucion en este espacio -&gt; vamos por metrica mixt" ID="ID_1061153029" CREATED="1456777736264" MODIFIED="1456778218924"/>
</node>
<node TEXT="Coherencia entre metrica transcripcional y GO (4.4)" POSITION="right" ID="ID_520927401" CREATED="1456778257120" MODIFIED="1456779578358">
<node TEXT="Interacting densities" ID="ID_858335920" CREATED="1456778344800" MODIFIED="1456778359014"/>
<node TEXT="KTA, zKTA" ID="ID_1185537005" CREATED="1456778360120" MODIFIED="1456778564453">
<node TEXT="figura geneX vs ontologias bpa, bpb y cc por tratamiento" ID="ID_1083558914" CREATED="1456779400504" MODIFIED="1456779424263"/>
<node TEXT="ariel: porque KTA control random PIN vs GO es TAN alto???? - descartar por ahora datos redes" ID="ID_754127809" CREATED="1456779424760" MODIFIED="1456779686935"/>
</node>
<node TEXT="fisher test (ariel)" ID="ID_519875816" CREATED="1456778564952" MODIFIED="1456779702294"/>
</node>
<node TEXT="Buscando Homogeneidad Biologica&#xa;( heterogeneidades en escala de clusters razonablemente extensos (DS1) )" POSITION="right" ID="ID_1743529587" CREATED="1456601397364" MODIFIED="1456601539871">
<edge COLOR="#ff00ff"/>
<node TEXT="concepto de lKTA&#xa;KTA en vecindades transcripcionales&#xa;para el analisis de TODOS los genes que se movieron en un tratamiento" ID="ID_1194508946" CREATED="1456601738316" MODIFIED="1456779982615">
<node TEXT="redes topologia&#xa;KMNN, (KNN, rcutoff(?): OUT)&#xa;Argumentar algo del valor 30 (Ariel)" ID="ID_1410452734" CREATED="1456605453502" MODIFIED="1456780177061">
<node TEXT="+fig1 lKTA30.R&#xa; + distribucion de grado&#xa; + betweennes" ID="ID_1520971845" CREATED="1456606058137" MODIFIED="1456606062324"/>
<node TEXT="Fig 2 lKTA30.R&#xa;+ nx vs ny (vecindad y vecinos anotados)" ID="ID_77968201" CREATED="1456606100212" MODIFIED="1456685062323"/>
<node TEXT="+wnyanotados vs wny: lo que le pasa&#xa;a un edge no predice comportamiento de vecindad" ID="ID_1269147884" CREATED="1456606063469" MODIFIED="1456686014395"/>
<node TEXT="lKTAanot vs wnyanot: relacion lineal" ID="ID_1801574071" CREATED="1456606214598" MODIFIED="1456606237819"/>
<node TEXT="lktaAnot vs ny: leve tendencia negativa" ID="ID_1196975713" CREATED="1456606238392" MODIFIED="1456606293716"/>
</node>
</node>
<node TEXT="metrica mixta.&#xa;idea: sobre la heterogeneidad topologica (transcripcional)&#xa;montar desorden de pesos a partir de info biologica" ID="ID_478424208" CREATED="1456601606646" MODIFIED="1456606997041">
<node TEXT="idea: hay muchas manera de mezclar metricas" ID="ID_1156613054" CREATED="1456605844988" MODIFIED="1456607941361">
<node TEXT="en el pasado el grupo exploro: mezcla convexa de distancias&#xa;dij = sqrt(alpha dijX ^2+ (1-alpha) dijGO^2)&#xa;busca el consenso: parametro global que parametriza de manera continua una metrica vs la otra [tesis_de_ariel, congresos_internacionales_poster_de_metrica_mixta]" ID="ID_1348779716" CREATED="1456607961312" MODIFIED="1456608292476"/>
<node TEXT="elegimos ahora no buscar el consenso sino penalizacion de correlaciones no soportadas por GOBP&#xa;wij &lt;-wij^stress &lt;- cor_ij ^(beta*stress)" ID="ID_28278735" CREATED="1456607904412" MODIFIED="1456781461680"/>
</node>
<node TEXT="redes pesos" ID="ID_118425612" CREATED="1456605555620" MODIFIED="1456605566630">
<node TEXT="&#xa;stress= lKTA_bkgd/lKTA_ij  (figura con distribucion de stress)&#xa;&#xa;wij &lt;- cor_ij ^ stress&#xa; queremos detectar heterogeneidades -&gt; transformamos los wij para que los strength de los nodos de la red tengan una distribucion de cola pesada eligiendo un beta tal que &#xa;&#xa;wij&apos; = wij^beta sea power law (por tratamiento)  (a la WGCNA)" ID="ID_1678676041" CREATED="1456605567832" MODIFIED="1456781164064">
<icon BUILTIN="help"/>
<node TEXT="para todos los genes del tratamiento" ID="ID_842369789" CREATED="1456780374728" MODIFIED="1456781120806"/>
<node TEXT="para clusters" ID="ID_1666132676" CREATED="1456781121336" MODIFIED="1456781141542"/>
</node>
</node>
</node>
<node TEXT="BHI analisis y algoritmo de merging" ID="ID_1024627094" CREATED="1457027910658" MODIFIED="1457027926512">
<node TEXT="(a modo de control):   (DS1)^2 (insideX)" ID="ID_681329268" CREATED="1456601415443" MODIFIED="1456601513468"/>
<node TEXT="penalizedMetric" ID="ID_1108544312" CREATED="1456609146569" MODIFIED="1456609166738">
<node TEXT="DTC. Penalizando algunos pesos de la red (los &quot;importantes&quot; en el sentido&#xa;de que aparecen como edges en la KMNN o KNN) que tienen much stress" ID="ID_507981328" CREATED="1456609167977" MODIFIED="1456609299177"/>
<node TEXT="graph clustering" ID="ID_1542212421" CREATED="1456609186862" MODIFIED="1456609201920">
<node TEXT="infomap" ID="ID_1004781034" CREATED="1456609305408" MODIFIED="1456609307527"/>
<node TEXT="cnm" ID="ID_1047856373" CREATED="1456609308231" MODIFIED="1456609318968"/>
</node>
</node>
</node>
<node TEXT="interpretacion biologica (Ariel)" ID="ID_653590148" CREATED="1456781397584" MODIFIED="1456781411534"/>
</node>
</node>
</map>
