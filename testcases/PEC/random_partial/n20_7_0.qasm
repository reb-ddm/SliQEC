OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
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
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
ccx q[2], q[11], q[10];
s q[16];
t q[6];
t q[12];
s q[2];
cx q[11], q[0];
s q[19];
cx q[5], q[1];
s q[1];
s q[14];
s q[1];
ccx q[16], q[11], q[18];
h q[17];
ccx q[4], q[16], q[8];
h q[15];
t q[7];
t q[0];
ccx q[15], q[0], q[11];
t q[6];
s q[0];
ccx q[5], q[9], q[8];
s q[5];
cx q[19], q[14];
h q[12];
ccx q[5], q[10], q[14];
ccx q[7], q[10], q[12];
ccx q[7], q[1], q[12];
t q[11];
cx q[3], q[16];
t q[17];
t q[14];
ccx q[10], q[17], q[9];
t q[11];
s q[10];
h q[15];
cx q[11], q[16];
s q[5];
t q[17];
cx q[16], q[18];
cx q[10], q[6];
h q[1];
s q[14];
cx q[12], q[3];
h q[5];
t q[8];
h q[5];
h q[18];
h q[12];
cx q[17], q[6];
h q[1];
t q[17];
h q[2];
h q[6];
h q[19];
h q[19];
cx q[13], q[6];
s q[14];
ccx q[10], q[1], q[17];
s q[8];
h q[1];
cx q[0], q[1];
tdg q[0];
tdg q[1];
sdg q[0];
cx q[0], q[1];
tdg q[3];
tdg q[2];
cx q[2], q[3];
tdg q[5];
cx q[4], q[5];
tdg q[5];
tdg q[4];
cx q[7], q[6];
tdg q[6];
cx q[7], q[6];
sdg q[6];
cx q[8], q[9];
h q[18];
s q[10];
h q[11];
ccx q[19], q[18], q[13];
ccx q[12], q[15], q[10];
ccx q[11], q[13], q[19];
t q[16];
s q[18];
t q[16];
cx q[14], q[15];
