OPENQASM 2.0;
include "qelib1.inc";
qreg q[47];
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
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
cx q[20], q[3];
cx q[25], q[37];
t q[5];
s q[32];
t q[5];
cx q[8], q[27];
s q[30];
ccx q[1], q[14], q[25];
t q[14];
ccx q[22], q[16], q[34];
s q[21];
t q[23];
ccx q[32], q[5], q[39];
cx q[3], q[11];
t q[37];
ccx q[35], q[34], q[20];
s q[19];
ccx q[23], q[28], q[15];
s q[38];
h q[24];
s q[25];
t q[38];
s q[10];
ccx q[6], q[32], q[14];
t q[34];
s q[7];
h q[26];
t q[6];
cx q[32], q[5];
ccx q[11], q[29], q[10];
s q[9];
cx q[27], q[22];
s q[3];
s q[33];
s q[6];
h q[39];
ccx q[34], q[15], q[22];
ccx q[31], q[17], q[6];
ccx q[39], q[25], q[30];
s q[9];
s q[16];
t q[3];
s q[18];
ccx q[27], q[0], q[12];
t q[1];
cx q[24], q[17];
cx q[1], q[34];
s q[25];
cx q[7], q[36];
s q[11];
cx q[4], q[31];
cx q[19], q[38];
ccx q[31], q[19], q[11];
h q[3];
s q[1];
cx q[33], q[11];
ccx q[1], q[29], q[38];
ccx q[7], q[13], q[29];
cx q[10], q[24];
cx q[12], q[37];
t q[26];
h q[2];
h q[1];
ccx q[33], q[22], q[30];
ccx q[7], q[30], q[3];
s q[17];
cx q[18], q[12];
s q[29];
cx q[20], q[28];
ccx q[22], q[31], q[10];
ccx q[15], q[22], q[0];
h q[35];
h q[31];
cx q[5], q[4];
cx q[1], q[29];
ccx q[21], q[17], q[13];
cx q[7], q[9];
cx q[39], q[0];
cx q[1], q[6];
t q[36];
h q[24];
ccx q[13], q[18], q[5];
t q[32];
h q[32];
s q[18];
h q[29];
s q[0];
s q[21];
t q[12];
cx q[11], q[19];
ccx q[23], q[31], q[33];
ccx q[21], q[27], q[3];
cx q[12], q[23];
cx q[23], q[18];
t q[21];
t q[3];
cx q[28], q[30];
cx q[37], q[0];
t q[39];
cx q[37], q[13];
s q[14];
ccx q[26], q[32], q[39];
cx q[14], q[36];
h q[2];
cx q[30], q[18];
ccx q[30], q[1], q[37];
t q[30];
h q[25];
s q[17];
h q[12];
ccx q[2], q[37], q[21];
ccx q[14], q[10], q[39];
s q[20];
h q[38];
h q[11];
ccx q[24], q[27], q[29];
cx q[23], q[35];
t q[23];
h q[4];
t q[31];
sdg q[0];
tdg q[1];
cx q[0], q[1];
sdg q[0];
y q[2];
tdg q[3];
cx q[3], q[4];
sdg q[3];
sdg q[6];
cx q[5], q[6];
tdg q[6];
cx q[5], q[6];
tdg q[5];
sdg q[7];
x q[9];
x q[11];
cx q[13], q[12];
x q[12];
sdg q[12];
cx q[13], q[12];
x q[12];
tdg q[15];
tdg q[14];
cx q[15], q[14];
sdg q[14];
cx q[16], q[17];
tdg q[17];
sdg q[17];
cx q[16], q[17];
tdg q[16];
cx q[19], q[18];
sdg q[18];
x q[18];
cx q[19], q[18];
x q[18];
s q[31];
h q[20];
h q[24];
h q[20];
ccx q[28], q[27], q[20];
h q[22];
ccx q[39], q[38], q[27];
h q[30];
cx q[37], q[27];
t q[23];
h q[32];
ccx q[20], q[26], q[30];
cx q[32], q[37];
h q[33];
ccx q[28], q[37], q[21];
s q[37];
ccx q[36], q[34], q[38];
t q[26];
ccx q[33], q[38], q[27];
s q[23];
cx q[40], q[2];
cx q[41], q[39];
cx q[42], q[8];
cx q[43], q[28];
cx q[44], q[37];
cx q[45], q[33];
cx q[46], q[38];
