OPENQASM 2.0;
include "qelib1.inc";
qreg q[9];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
t q[1];
cx q[5], q[3];
t q[2];
ccx q[6], q[3], q[1];
cx q[2], q[1];
h q[5];
cx q[7], q[0];
t q[4];
ccx q[2], q[5], q[6];
t q[0];
s q[4];
t q[3];
h q[7];
h q[7];
t q[2];
t q[1];
h q[3];
s q[5];
ccx q[2], q[7], q[6];
s q[4];
s q[0];
cx q[4], q[2];
cx q[4], q[0];
h q[0];
cx q[0], q[1];
sdg q[0];
x q[2];
cx q[2], q[3];
sdg q[3];
cx q[2], q[3];
x q[2];
s q[6];
s q[5];
ccx q[7], q[4], q[5];
ccx q[6], q[7], q[4];
cx q[8], q[4];
