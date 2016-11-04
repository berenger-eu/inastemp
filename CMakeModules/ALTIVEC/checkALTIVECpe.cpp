
#include <altivec.h>

int main() {
	__vector double tx;
    __vector double ty;
	tx = ty + tx;
	tx = ty * tx;
	tx = ty - tx;
	tx = ty / tx;
	tx += ty;
	tx -= ty;
	tx *= ty;
	tx /= ty;
  return 0;
}
