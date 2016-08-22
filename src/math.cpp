#include <node.h>

static void init(v8::Handle<v8::Object> exports, v8::Handle<v8::Object> module)
{
	//nothing to register, only for require('math')
}

NODE_MODULE(math, init)