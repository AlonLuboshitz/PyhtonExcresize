class A:
    pass
class B(A):
    pass
    class C(A):
        def __init__(self) -> None:
            pass
c = B.C()
print(type(c))
