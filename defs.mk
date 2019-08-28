INSTALL=/usr/bin/install

PREFIX=/usr/local


PYTHON2.7=$(shell if [ -x /bin/python2.7 ] ; then echo "/bin/python2.7" ; elif [ -x /usr/bin/python2.7 ] ; then echo "/usr/bin/python2.7" ; elif [ -x /usr/local/bin/python2.7 ] ; then echo "/usr/local/bin/python2.7" ; elif [ -x /opt/bin/python2.7 ] ; then echo "/opt/bin/python2.7" ; else echo "/none" ; fi )

PYTHON2.6=$(shell if [ -x /bin/python2.6 ] ; then echo "/bin/python2.6" ; elif [ -x /usr/bin/python2.6 ] ; then echo "/usr/bin/python2.6" ; elif [ -x /usr/local/bin/python2.6 ] ; then echo "/usr/local/bin/python2.6" ; elif [ -x /opt/bin/python2.6 ] ; then echo "/opt/bin/python2.6" ; else echo "/none" ; fi )

PYTHON3.4=$(shell if [ -x /bin/python3.4 ] ; then echo "/bin/python3.4" ; elif [ -x /usr/bin/python3.4 ] ; then echo "/usr/bin/python3.4" ; elif [ -x /usr/local/bin/python3.4 ] ; then echo "/usr/local/bin/python3.4" ; elif [ -x /opt/bin/python3.4 ] ; then echo "/opt/bin/python3.4" ; else echo "/none" ; fi )


PYTHON3.6=$(shell if [ -x /bin/python3.6 ] ; then echo "/bin/python3.6" ; elif [ -x /usr/bin/python3.6 ] ; then echo "/usr/bin/python3.6" ; elif [ -x /usr/local/bin/python3.6 ] ; then echo "/usr/local/bin/python3.6" ; elif [ -x /opt/bin/python3.6 ] ; then echo "/opt/bin/python3.6" ; elif [ -x /bin/python36 ] ; then echo "/bin/python36" ; elif [ -x /usr/bin/python36 ] ; then echo "/usr/bin/python36" ; elif [ -x /usr/local/bin/python36 ] ; then echo "/usr/local/bin/python36" ; elif [ -x /opt/bin/python36 ] ; then echo "/opt/bin/python36" ; elif [ -x /opt/rh/rh-python36/root/bin/python3.6 ] ; then echo "/opt/rh/rh-python36/root/bin/python3.6" ;else echo "/none" ; fi )

PYTHON3.7=$(shell if [ -x /bin/python3.7 ] ; then echo "/bin/python3.7" ; elif [ -x /usr/bin/python3.7 ] ; then echo "/usr/bin/python3.7" ; elif [ -x /usr/local/bin/python3.7 ] ; then echo "/usr/local/bin/python3.7" ; elif [ -x /opt/bin/python3.7 ] ; then echo "/opt/bin/python3.7" ; elif [ -x /bin/python37 ] ; then echo "/bin/python37" ; elif [ -x /usr/bin/python37 ] ; then echo "/usr/bin/python37" ; elif [ -x /usr/local/bin/python37 ] ; then echo "/usr/local/bin/python37" ; elif [ -x /opt/bin/python37 ] ; then echo "/opt/bin/python37" ; elif [ -x /opt/rh/rh-python37/root/bin/python3.7 ] ; then echo "/opt/rh/rh-python37/root/bin/python3.7" ;else echo "/none" ; fi )

PYTHON3.8=$(shell if [ -x /bin/python3.8 ] ; then echo "/bin/python3.8" ; elif [ -x /usr/bin/python3.8 ] ; then echo "/usr/bin/python3.8" ; elif [ -x /usr/local/bin/python3.8 ] ; then echo "/usr/local/bin/python3.8" ; elif [ -x /opt/bin/python3.8 ] ; then echo "/opt/bin/python3.8" ; elif [ -x /bin/python38 ] ; then echo "/bin/python38" ; elif [ -x /usr/bin/python38 ] ; then echo "/usr/bin/python38" ; elif [ -x /usr/local/bin/python38 ] ; then echo "/usr/local/bin/python38" ; elif [ -x /opt/bin/python38 ] ; then echo "/opt/bin/python38" ; elif [ -x /opt/rh/rh-python38/root/bin/python3.8 ] ; then echo "/opt/rh/rh-python38/root/bin/python3.8" ;else echo "/none" ; fi )



DEFAULTPY=$(shell if [ -x /bin/python ] ; then echo "/bin/python" ; elif [ -x /usr/bin/python ] ; then echo "/usr/bin/python" ; elif [ -x /usr/local/bin/python ] ; then echo "/usr/local/bin/python" ; elif [ -x /opt/bin/python ] ; then echo "/opt/bin/python" ; else echo "/none" ; fi )

