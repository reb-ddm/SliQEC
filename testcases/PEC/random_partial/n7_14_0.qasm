OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
cx q[2], q[3];
s q[6];
cx q[2], q[0];
t q[5];
ccx q[1], q[4], q[6];
ccx q[2], q[5], q[0];
s q[5];
h q[6];
ccx q[4], q[5], q[6];
cx q[6], q[4];
cx q[2], q[1];
h q[6];
cx q[2], q[6];
cx q[0], q[5];
h q[0];
s q[3];
h q[0];
cx q[3], q[4];
s q[0];
t q[0];
s q[3];
sdg q[0];
y q[2];
cx q[4], q[3];
s q[4];
ccx q[3], q[6], q[4];
cx q[5], q[3];
