#pragma once

const char GREEN[] = "\033[1;32m";
const char RED[] = "\033[1;31m";
const char WHITE[] = "\033[1;37m";
const char NORM[] = "\033[0m";
const float EPSILON = 0.000001;

#define TEST_ASSERT(what,op,ref)\
{\
	auto result = (what);\
	auto comp = (ref);\
	\
	if(!(result op comp))\
    {\
        std::cout << RED << "Assertion failed:" << std::endl;\
		std::cout << #what << " = " << result << std::endl;\
		std::cout << #ref << " = " << comp << NORM << std::endl;\
        assert(what op ref);\
    }\
}