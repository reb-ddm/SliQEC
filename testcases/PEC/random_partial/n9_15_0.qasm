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
h q[8];
s q[6];
ccx q[1], q[3], q[6];
cx q[7], q[4];
t q[8];
h q[1];
cx q[6], q[8];
s q[2];
t q[1];
h q[4];
s q[5];
h q[2];
s q[8];
cx q[0], q[4];
h q[1];
ccx q[6], q[1], q[3];
ccx q[8], q[0], q[6];
ccx q[1], q[8], q[7];
h q[6];
ccx q[0], q[6], q[2];
ccx q[6], q[0], q[3];
ccx q[4], q[5], q[3];
h q[4];
h q[4];
t q[3];
s q[4];
s q[8];
h q[3];
tdg q[0];
tdg q[1];
cx q[1], q[0];
sdg q[0];
sdg q[2];
ccx q[5], q[8], q[4];
s q[7];
ccx q[8], q[5], q[7];
cx q[5], q[8];
h q[6];
