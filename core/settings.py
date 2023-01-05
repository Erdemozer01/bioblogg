import os
from pathlib import Path
from django.shortcuts import redirect
from django.urls import reverse_lazy

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/4.1/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'django-insecure-_*6$iz2emf&x(*^w5#)^j7o7cyx^@d%z+2_@%e0(#(qh786$*_'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

ALLOWED_HOSTS = ['bioblog.pythonanywhere.com', 'localhost', '127.0.0.1']

# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'blog.apps.BlogConfig',
    'accounts.apps.AccountsConfig',
    'django_plotly_dash.apps.DjangoPlotlyDashConfig',
    'dashboard.apps.DashboardConfig',
    'django_bootstrap5',
    'bootstrap4',
    'django_cleanup',
    'hitcount',
    'ckeditor',
    'ckeditor_uploader',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'core.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [BASE_DIR / 'templates']
        ,
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'core.wsgi.application'

# Database
# https://docs.djangoproject.com/en/4.1/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / 'db.sqlite3'
    }
}

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

# Internationalization
# https://docs.djangoproject.com/en/4.1/topics/i18n/

LANGUAGE_CODE = 'tr'

TIME_ZONE = 'Europe/Istanbul'

USE_I18N = True

USE_TZ = True

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/4.1/howto/static-files/

STATIC_URL = 'static/'

STATICFILES_DIRS = [
    BASE_DIR / "static"
]

MEDIA_ROOT = os.path.join(BASE_DIR, 'media')

MEDIA_URL = "/media/"

LOGIN_REDIRECT_URL = reverse_lazy('blog:anasayfa')

LOGOUT_REDIRECT_URL = reverse_lazy('blog:anasayfa')

DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'

CKEDITOR_UPLOAD_PATH = "uploads/"

X_FRAME_OPTIONS = 'SAMEORIGIN'

CKEDITOR_CONFIGS = {
    'default': {
        'toolbar': 'full',
        'height': 'auto',
        'width': 'auto',
    },
    'blog': {
        'height': 'auto',
        'width': 'auto',
        'toolbar_Basic': [
            ['Source', '-', 'Bold', 'Italic']
        ],
        'toolbar_BlogToolbarConfig': [
            {'name': 'styles', 'items': ['Format', 'Font', 'FontSize']},
            {'name': 'basicstyles',
             'items': ['JustifyLeft', 'JustifyCenter', 'JustifyRight', 'JustifyBlock', '-', 'Bold', 'Italic',
                       'Underline', 'Strike', 'Subscript', 'Superscript', 'RemoveFormat', 'TextColor', 'Maximize']},
            {'name': 'paragraph',
             'items': ['Blockquote', '-', 'NumberedList', 'BulletedList', 'ExportPdf', ]},
            {'name': 'insert', 'items': ['Flash', 'Table', 'Smiley', 'SpecialChar']},
            {'name': 'yourcustomtools', 'items': ['Undo', 'Redo']},
            {'name': 'links', 'items': ['Link', 'Unlink', '-', ]},
        ],
        'toolbar': 'BlogToolbarConfig',
        'tabSpaces': 4,
        'mathJaxLib': '//cdnjs.cloudflare.com/ajax/libs/mathjax/3.0.0/es5/latest?tex-mml-chtml.js',
    },
}
