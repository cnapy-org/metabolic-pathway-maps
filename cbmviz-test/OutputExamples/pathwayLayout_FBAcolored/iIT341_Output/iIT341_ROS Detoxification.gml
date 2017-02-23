graph [
directed 1
version 2
node [
id 0
label "M_h2o2_c0"
graphics [
x -0.75
y 2.75
z 0
h 2.5
w 2.5
d 1
type "rectangle"
width 0.12
fill "#67a9cf"
outline "#000000"
     ]
]
node [
id 1
label "R_SPODM1"
graphics [
x 2.25
y 7.5
z 0
h 1
w 1
d 1
type "rectangle"
width 0.12
fill "#fc9272"
outline "#000000"
     ]
]
node [
id 2
label "M_o2s_c2"
graphics [
x 2.25
y 12.25
z 0
h 2.5
w 2.5
d 1
type "rectangle"
width 0.12
fill "#67a9cf"
outline "#000000"
     ]
]
node [
id 3
label "M_ficytcc553_c3"
graphics [
x -6.75
y -6.75
z 0
h 2.5
w 2.5
d 1
type "rectangle"
width 0.12
fill "#67a9cf"
outline "#000000"
     ]
]
node [
id 4
label "M_focytcc553_c4"
graphics [
x -6.75
y 2.75
z 0
h 2.5
w 2.5
d 1
type "rectangle"
width 0.12
fill "#67a9cf"
outline "#000000"
     ]
]
node [
id 5
label "M_o2_c5"
graphics [
x 2.25
y -6.75
z 0
h 2.5
w 2.5
d 1
type "rectangle"
width 0.12
fill "#67a9cf"
outline "#000000"
     ]
]
node [
id 6
label "R_CAT6"
graphics [
x -0.75
y -2
z 0
h 1
w 1
d 1
type "rectangle"
width 0.12
fill "#fc9272"
outline "#000000"
     ]
]
node [
id 7
label "R_CCP7"
graphics [
x -6.75
y -2
z 0
h 1
w 1
d 1
type "rectangle"
width 0.12
fill "#fc9272"
outline "#000000"
     ]
]
node [
id 8
label "ROS Detoxification8"
graphics [
x 0
y 0
z 0
h 5
w 5
d 1
type "rectangle"
width 0.12
fill "#ebebeb"
outline "#000000"
     ]
]
edge [
source 2
target 1
id 0
label ""
graphics [
type "line"
arrow "last"
width 0.1
Line [
]
]
]
edge [
source 1
target 5
id 1
label ""
graphics [
type "line"
arrow "last"
width 0.1
Line [
point [
x 2.25
y 7.5
z 0
]
point [
x 5.25
y 2.75
z 0
]
point [
x 5.25
y -2
z 0
]
point [
x 2.25
y -6.75
z 0
]
]
]
]
edge [
source 6
target 5
id 2
label ""
graphics [
type "line"
arrow "last"
width 0.1
Line [
]
]
]
edge [
source 1
target 0
id 3
label ""
graphics [
type "line"
arrow "last"
width 0.1
Line [
]
]
]
edge [
source 7
target 3
id 4
label ""
graphics [
type "line"
arrow "last"
width 0.1
Line [
]
]
]
edge [
source 4
target 7
id 5
label ""
graphics [
type "line"
arrow "last"
width 0.1
Line [
]
]
]
edge [
source 0
target 6
id 6
label ""
graphics [
type "line"
arrow "last"
width 0.1
Line [
]
]
]
edge [
source 0
target 7
id 7
label ""
graphics [
type "line"
arrow "last"
width 0.1
Line [
]
]
]
]
