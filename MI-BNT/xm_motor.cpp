#include "xm_motor.h"
#include <Arduino.h>
MCP2515 mcp2515(10);
struct can_frame canMsg;
uint32_t ExtId;             //定义can扩展id
uint8_t rx_data[8];         //接收数据
uint32_t Motor_Can_ID;      //接收数据电机ID
static uint8_t byte_ls[4];  //转换临时数据
uint8_t tx_data[8];         //can写入的数据

int filter(int queue[], char n)  //数值过滤
{
  int sum = 0;
  byte i;
  int maxsz = queue[0];  //寻找最大值
  int minsz = queue[0];  //寻找最小值
  for (i = 0; i < n; i++) {
    if (maxsz < queue[i]) {
      maxsz = queue[i];  //寻找最大值和最小值
    }
    if (minsz > queue[i]) {
      minsz = queue[i];
    }
  }

  for (i = 0; i < n; i++) {
    sum += queue[i];
  }
  sum = sum - maxsz - minsz;  //去除最大值和最小值
  return (sum / (n - 2));
}  //数值过滤

int check(byte ao_port, byte n)  //采样
{

  int check_date[n];  //定义采样数组
  for (byte i = 0; i < n; ++i) {
    check_date[i] = analogRead(ao_port);  //获得指定传感器数据
  }

  int vvvv = filter(check_date, n);
  return vvvv;

}  //采样

void xm_can_start() {

  mcp2515.reset();
  mcp2515.setBitrate(CAN_1000KBPS, MCP_8MHZ);
  mcp2515.setNormalMode();
  delay(4);
}
static uint8_t* Float_to_Byte(float f)  //float分解成四个byte数据
{
  unsigned long longdata = 0;
  longdata = *(unsigned long*)&f;
  byte_ls[0] = (longdata & 0xFF000000) >> 24;
  byte_ls[1] = (longdata & 0x00FF0000) >> 16;
  byte_ls[2] = (longdata & 0x0000FF00) >> 8;
  byte_ls[3] = (longdata & 0x000000FF);
  return byte_ls;
}
static float uint16_to_float(uint16_t x, float x_min, float x_max, int bits)  //把uint 16位数据变成浮点数 用在接受数据的处理上
{
  uint32_t span = (1 << bits) - 1;
  float offset = x_max - x_min;
  return offset * x / span + x_min;
}
static int float_to_uint(float x, float x_min, float x_max, int bits)  //把浮点数转换成uint_16 用在位置 扭矩 上面
{
  float span = x_max - x_min;
  float offset = x_min;
  if (x > x_max) x = x_max;
  else if (x < x_min) x = x_min;
  return (int)((x - offset) * ((float)((1 << bits) - 1)) / span);
}
void exid_count(uint8_t Communication_Type, uint16_t msid, uint8_t can_id)  //计算扩展ExtId,Communication_Type通信类型，msid主canid，
{
  uint8_t msid_l = msid;
  uint8_t msid_h = msid >> 8;
  uint32_t di_data = ((0xFFFFFFFF & Communication_Type) << 24) | 0x00FFFFFF;  //求出高32位
  uint32_t di_datab = ((0xFFFFFFFF & msid_h) << 16) | 0xFF00FFFF;             //求出高32位
  uint32_t di_datac = ((0xFFFFFFFF & msid_l) << 8) | 0xFFFF00FF;              //求出高32位
  uint32_t di_datad = (0xFFFFFFFF & can_id) | 0xFFFFFF00;                     //求出高32位
  ExtId = (di_data & di_datab & di_datac & di_datad);
}  //计算扩展ExtId,Communication_Type通信类型，msid主canid，

void data_count_dcs(uint16_t Index, float Value, char Value_type) {
  //计算can在 单参数写入，通信类型 12下发送的8位数据，Index 是命令类型0: 运控模式1: 位置模式2: 速度模式3: 电流模式Value是0 值，Value_type是数据类型，浮点数用f非浮点用s
  //速度数值要 注明浮点数 f ,
  //写入扭矩 n


  canMsg.data[0] = Index;
  canMsg.data[1] = Index >> 8;
  canMsg.data[2] = 0x00;
  canMsg.data[3] = 0x00;
  if (Value_type == 'f') {
    Float_to_Byte(Value);
    canMsg.data[4] = byte_ls[3];
    canMsg.data[5] = byte_ls[2];
    canMsg.data[6] = byte_ls[1];
    canMsg.data[7] = byte_ls[0];
  } else if (Value_type == 's') {
    canMsg.data[4] = (uint8_t)Value;
    canMsg.data[5] = 0x00;
    canMsg.data[6] = 0x00;
    canMsg.data[7] = 0x00;
  }
}  //计算can在 单参数写入，通信类型 12下发送的8位数据，Index 是命令类型Value是值，Value_type是数据类型，浮点数用f非浮点用s
void data_count_zero()  //can数据置零
{
  canMsg.data[0] = 0x00;
  canMsg.data[1] = 0x00;
  canMsg.data[2] = 0x00;
  canMsg.data[3] = 0x00;
  canMsg.data[4] = 0x00;
  canMsg.data[5] = 0x00;
  canMsg.data[6] = 0x00;
  canMsg.data[7] = 0x00;
}  //can数据置零

void motor_enable(uint8_t id)  //电机使能 电机canid
{
  exid_count(3, Master_CAN_ID, id);
  canMsg.can_id = ExtId | CAN_EFF_FLAG;
  canMsg.can_dlc = 8;
  data_count_zero();
  mcp2515.sendMessage(MCP2515::TXB1, &canMsg);
  delay(4);
}  //电机使能 电机canid
void motor_mode(uint8_t id, char type)  //电机运行模式 电机canid  模式值1位置模式2速度模式 3 电流模式0运控模式
{
  exid_count(0x12, Master_CAN_ID, id);
  canMsg.can_id = ExtId | CAN_EFF_FLAG;  //ExtId  0x12000001
  canMsg.can_dlc = 8;
  data_count_dcs(0x7005, type, 's');
  mcp2515.sendMessage(MCP2515::TXB1, &canMsg);
  delay(4);
}

void motor_speed_cmd(uint8_t id, float speed_ref) {  //设置速度模式下的参数转速
  exid_count(0x12, Master_CAN_ID, id);
  canMsg.can_id = ExtId | CAN_EFF_FLAG;  //ExtId  0x12000001
  canMsg.can_dlc = 8;
  data_count_dcs(0x700A, speed_ref, 'f');
  mcp2515.sendMessage(MCP2515::TXB1, &canMsg);
//  delay(1);
}

void motor_pos_zero(uint8_t id)  //位置回0
{
  exid_count(6, Master_CAN_ID, id);
  canMsg.can_id = ExtId | CAN_EFF_FLAG;
  canMsg.can_dlc = 8;
  data_count_zero();
  canMsg.data[0] = 1;
  mcp2515.sendMessage(MCP2515::TXB1, &canMsg);
  delay(4);
}  //位置回0

void motor_set_zero(uint8_t id)  //位置置0
{
  exid_count(6, Master_CAN_ID, id);
  canMsg.can_id = ExtId | CAN_EFF_FLAG;
  canMsg.can_dlc = 8;
  data_count_zero();
  canMsg.data[0] = 0;
  mcp2515.sendMessage(MCP2515::TXB1, &canMsg);
  delay(4);
}  //位置置0

void motor_pos_cmd(uint8_t id, float pos) {  //设置位置模式下的位置
  exid_count(0x12, Master_CAN_ID, id);
  canMsg.can_id = ExtId | CAN_EFF_FLAG;  //ExtId  0x12000001
  canMsg.can_dlc = 8;
  data_count_dcs(0x7016, pos, 'f');
  mcp2515.sendMessage(MCP2515::TXB1, &canMsg);
  //  delay(1);
//  delayMicroseconds(5);
}

void motor_pow_cmd(uint8_t id, float lim_spd, float lim_cur, float kp = 0.125, float ki = 0.0158) {  //设置速度 电流限制，kp,kd
  exid_count(0x12, Master_CAN_ID, id);
  canMsg.can_id = ExtId | CAN_EFF_FLAG;  //ExtId  0x12000001
  canMsg.can_dlc = 8;
  data_count_dcs(0x7017, lim_spd, 'f');
  mcp2515.sendMessage(MCP2515::TXB1, &canMsg);
  delay(4);
  data_count_dcs(0x7018, lim_cur, 'f');
  mcp2515.sendMessage(MCP2515::TXB1, &canMsg);
  delay(4);
  data_count_dcs(0x7010, kp, 'f');
  mcp2515.sendMessage(MCP2515::TXB1, &canMsg);
  delay(4);
  data_count_dcs(0x7011, ki, 'f');
  mcp2515.sendMessage(MCP2515::TXB1, &canMsg);
  delay(4);
}

//int 
