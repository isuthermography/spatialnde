
include defs.mk

SUBDIRS= dataguzzler

all:
	@for i in $(SUBDIRS) ; do if [ -d $$i ] && [ -f $$i/Makefile ] ; then $(MAKE) $(MFLAGS) -C $$i ; fi done
ifneq ($(PYTHON2.6), /none)
	$(PYTHON2.6) ./setup.py build
endif
ifneq ($(PYTHON2.7), /none)
	$(PYTHON2.7) ./setup.py build 
endif
ifneq ($(PYTHON3.4), /none)
	$(PYTHON3.4) ./setup.py build
endif
ifneq ($(PYTHON3.6), /none)
	$(PYTHON3.6) ./setup.py build 
endif
ifneq ($(PYTHON3.7), /none)
	$(PYTHON3.7) ./setup.py build 
endif
ifneq ($(PYTHON3.8), /none)
	$(PYTHON3.8) ./setup.py build 
endif
ifneq ($(DEFAULTPY), /none)
	$(DEFAULTPY) ./setup.py build
endif

clean:
	@for i in $(SUBDIRS) ; do if [ -d $$i ] && [ -f $$i/Makefile ] ; then $(MAKE) $(MFLAGS) -C $$i clean ; fi done
	#python setup.py clean
	rm -rf build/
	rm -f *.pyc *~ */*~ */*.pyc *.bak */*.so */*.dll */*.bak *.aux *.log *.pdf spatial/*~ demos/*~ demos/*.pyc
	rm -f */*/*~ */*/*.pyc */*/*.bak


install:
	@for i in $(SUBDIRS) ; do if [ -d $$i ] && [ -f $$i/Makefile ] ; then $(MAKE) $(MFLAGS) -C $$i install ; fi done
ifneq ($(PYTHON2.6), /none)
	$(PYTHON2.6) ./setup.py install --prefix=$(PREFIX) 
endif
ifneq ($(PYTHON2.7), /none)
	$(PYTHON2.7) ./setup.py install --prefix=$(PREFIX) 
endif
ifneq ($(PYTHON3.4), /none)
	$(PYTHON3.4) ./setup.py install --prefix=$(PREFIX) 
endif
ifneq ($(PYTHON3.6), /none)
	$(PYTHON3.6) ./setup.py install --prefix=$(PREFIX) 
endif
ifneq ($(PYTHON3.7), /none)
	$(PYTHON3.7) ./setup.py install --prefix=$(PREFIX) 
endif
ifneq ($(PYTHON3.8), /none)
	$(PYTHON3.8) ./setup.py install --prefix=$(PREFIX) 
endif
ifneq ($(DEFAULTPY), /none)
	$(DEFAULTPY) ./setup.py install --prefix=$(PREFIX) 
	#$(DEFAULTPY) ./setup.py install_data # --prefix=$(PREFIX) 
endif
