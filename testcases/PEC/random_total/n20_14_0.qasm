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
cx q[10], q[19];
t q[10];
h q[6];
t q[5];
s q[17];
s q[15];
s q[2];
ccx q[8], q[3], q[1];
cx q[2], q[12];
cx q[15], q[4];
s q[1];
h q[4];
h q[12];
ccx q[3], q[9], q[2];
ccx q[18], q[3], q[6];
t q[12];
s q[16];
cx q[8], q[10];
s q[12];
t q[5];
s q[2];
s q[16];
t q[14];
h q[10];
ccx q[13], q[9], q[14];
h q[4];
t q[5];
cx q[7], q[19];
h q[4];
cx q[13], q[3];
t q[7];
h q[10];
s q[6];
s q[19];
ccx q[10], q[15], q[19];
ccx q[15], q[16], q[8];
t q[0];
ccx q[19], q[5], q[8];
s q[16];
ccx q[3], q[9], q[14];
s q[3];
h q[14];
ccx q[16], q[17], q[18];
cx q[13], q[10];
ccx q[10], q[7], q[19];
cx q[4], q[16];
h q[15];
s q[0];
h q[17];
h q[4];
h q[1];
t q[9];
s q[2];
h q[7];
s q[11];
cx q[16], q[11];
cx q[6], q[12];
h q[19];
s q[18];
h q[14];