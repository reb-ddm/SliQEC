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
cx q[5], q[19];
t q[4];
t q[7];
s q[9];
s q[2];
ccx q[7], q[1], q[10];
ccx q[3], q[10], q[16];
ccx q[8], q[17], q[0];
ccx q[13], q[18], q[8];
h q[14];
cx q[14], q[10];
t q[4];
cx q[3], q[9];
ccx q[4], q[15], q[2];
ccx q[13], q[11], q[0];
h q[9];
h q[18];
h q[3];
cx q[8], q[11];
t q[7];
t q[8];
h q[3];
t q[8];
cx q[1], q[9];
t q[13];
s q[8];
cx q[16], q[17];
t q[0];
t q[19];
cx q[9], q[8];
h q[7];
t q[10];
h q[12];
cx q[6], q[0];
cx q[8], q[16];
ccx q[16], q[11], q[6];
ccx q[7], q[18], q[17];
ccx q[2], q[3], q[18];
h q[19];
t q[5];
cx q[17], q[0];
s q[4];
cx q[16], q[12];
ccx q[1], q[7], q[5];
ccx q[5], q[1], q[17];
ccx q[10], q[11], q[0];
s q[2];
ccx q[10], q[5], q[18];
h q[2];
ccx q[10], q[14], q[0];
ccx q[19], q[9], q[16];
cx q[5], q[16];
t q[11];
cx q[17], q[16];
cx q[14], q[18];
cx q[2], q[18];
t q[15];
cx q[1], q[13];
s q[15];
ccx q[7], q[0], q[17];
sdg q[0];
tdg q[2];
sdg q[4];
tdg q[5];
cx q[5], q[4];
sdg q[4];
sdg q[6];
cx q[6], q[7];
x q[8];
s q[11];
ccx q[13], q[11], q[16];
cx q[11], q[16];
h q[14];
ccx q[17], q[13], q[16];
cx q[14], q[18];
s q[19];
h q[18];
cx q[18], q[19];
ccx q[17], q[13], q[16];
