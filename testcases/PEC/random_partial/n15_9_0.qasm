OPENQASM 2.0;
include "qelib1.inc";
qreg q[15];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
cx q[4], q[7];
ccx q[12], q[6], q[3];
t q[3];
cx q[1], q[2];
h q[13];
cx q[5], q[1];
h q[11];
h q[13];
s q[12];
cx q[4], q[6];
s q[8];
cx q[1], q[14];
ccx q[2], q[1], q[5];
s q[9];
t q[14];
t q[3];
h q[6];
s q[12];
t q[7];
s q[6];
t q[1];
ccx q[9], q[4], q[2];
s q[13];
s q[12];
cx q[0], q[6];
t q[9];
h q[6];
t q[8];
ccx q[12], q[0], q[4];
s q[0];
ccx q[0], q[8], q[12];
t q[4];
ccx q[7], q[14], q[0];
h q[8];
t q[12];
ccx q[9], q[14], q[0];
t q[7];
cx q[2], q[13];
s q[0];
cx q[14], q[13];
h q[14];
t q[10];
cx q[14], q[2];
h q[1];
s q[13];
cx q[1], q[0];
tdg q[0];
tdg q[3];
tdg q[2];
cx q[3], q[2];
tdg q[2];
y q[4];
sdg q[5];
s q[9];
ccx q[11], q[8], q[14];
s q[11];
s q[7];
h q[13];
h q[12];
s q[10];
h q[7];
