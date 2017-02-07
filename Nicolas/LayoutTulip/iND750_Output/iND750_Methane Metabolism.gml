graph [
directed 1
version 2
node [
id 0
label "M_h2o_c0"
graphics [
x 1.75
y 12.25
z 0
h 1.5
w 1.5
d 1
type "rectangle"
width 0.12
fill "#ffff99"
outline "#000000"
     ]
]
node [
id 1
label "M_co2_c1"
graphics [
x -1.25
y -6.75
z 0
h 1.5
w 1.5
d 1
type "rectangle"
width 0.12
fill "#ffff99"
outline "#000000"
     ]
]
node [
id 2
label "M_nad_c2"
graphics [
x -7.25
y 2.75
z 0
h 1.5
w 1.5
d 1
type "rectangle"
width 0.12
fill "#ffff99"
outline "#000000"
     ]
]
node [
id 3
label "M_gthrd_c3"
graphics [
x -1.25
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
id 4
label "M_Sfglutth_c4"
graphics [
x 7.75
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
id 5
label "M_h_c5"
graphics [
x 4.75
y 2.75
z 0
h 1.5
w 1.5
d 1
type "rectangle"
width 0.12
fill "#ffff99"
outline "#000000"
     ]
]
node [
id 6
label "R_FDH6"
graphics [
x 1.75
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
label "R_SFGTHi7"
graphics [
x 4.75
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
id 8
label "M_nadh_c8"
graphics [
x 4.75
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
id 9
label "M_for_c9"
graphics [
x 10.75
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
id 10
label "Methane Metabolism10"
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
source 6
target 8
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
source 9
target 6
id 1
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
target 9
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
source 4
target 7
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
source 2
target 6
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
source 6
target 1
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
edge [
source 7
target 5
id 8
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
