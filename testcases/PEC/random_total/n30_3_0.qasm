OPENQASM 2.0;
include "qelib1.inc";
qreg q[30];
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
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
cx q[25], q[14];
s q[8];
ccx q[11], q[12], q[27];
h q[18];
cx q[14], q[21];
h q[5];
h q[25];
ccx q[25], q[8], q[0];
cx q[12], q[6];
h q[19];
s q[25];
h q[13];
cx q[10], q[28];
cx q[25], q[23];
cx q[1], q[8];
cx q[21], q[18];
ccx q[0], q[29], q[2];
h q[4];
ccx q[17], q[7], q[5];
t q[13];
s q[16];
t q[29];
cx q[7], q[20];
h q[8];
s q[19];
s q[7];
s q[0];
h q[14];
t q[9];
ccx q[0], q[12], q[14];
ccx q[24], q[20], q[5];
h q[5];
h q[16];
cx q[11], q[2];
t q[0];
cx q[23], q[16];
s q[8];
cx q[1], q[13];
cx q[15], q[5];
s q[20];
h q[27];
ccx q[16], q[3], q[6];
ccx q[29], q[12], q[10];
t q[7];
ccx q[1], q[10], q[21];
t q[28];
h q[28];
ccx q[18], q[26], q[25];
s q[3];
cx q[24], q[0];
cx q[11], q[8];
cx q[28], q[5];
h q[22];
ccx q[5], q[27], q[13];
s q[9];
ccx q[17], q[18], q[1];
h q[7];
h q[18];
t q[25];
t q[14];
ccx q[9], q[8], q[16];
t q[11];
s q[2];
cx q[17], q[3];
h q[4];
s q[10];
cx q[6], q[17];
t q[23];
t q[20];
t q[8];
ccx q[27], q[23], q[0];
ccx q[25], q[28], q[29];
h q[21];
cx q[21], q[29];
h q[18];
cx q[14], q[25];
h q[23];
h q[23];
cx q[26], q[20];
s q[12];
s q[13];
ccx q[1], q[11], q[13];
t q[9];
s q[2];
ccx q[14], q[17], q[29];
t q[29];
ccx q[6], q[17], q[21];
h q[26];
s q[14];
ccx q[13], q[4], q[17];
