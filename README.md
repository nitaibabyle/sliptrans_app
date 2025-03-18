This work is inspired by Lianqiang Han and Pierre Gergondet's work.
The IKQP is based on the Lianqiang Han's C++ work. I reconstucted the class of the solver, and achieves a generalized API.
The QP solver is based on Pierre's work, and Han gives a convenient port.
The robot2slip is achieved by forward kinematics, which is solved by pinocchio, a open source dynamics lib.
The slip2robot is achieved by inverse kinematics, which is solved by a QP problem. 
test.cpp shows a example of the usage for this transation tool. You can use the functions to refer to the code.
