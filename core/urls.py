from django.contrib import admin
from django.urls import path, include
from django.views import generic
from django.conf import settings
from django.conf.urls.static import static
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from blog.views import ArchiveView

urlpatterns = [
                  path('admin/', admin.site.urls),
                  path('accounts/', include('django.contrib.auth.urls')),
                  path('laboratuvar/bioinformatic/', include('django_plotly_dash.urls', namespace='bioinformatic_dash')),
                  path('laboratuvar/biostatistic/', include('django_plotly_dash.urls',  namespace='biostatistic_dash')),
                  path("ckeditor5/", include('django_ckeditor_5.urls')),
                  path('blog/', include("blog.urls", namespace='blog')),
                  path('accounts/', include("accounts.urls", namespace='accounts')),
                  path('dashboard/', include("dashboard.urls", namespace='dashboard')),
                  path('laboratuvar/bioinformatic/', include('bioinformatic.urls', namespace='bioinformatic')),
                  path('laboratuvar/biostatistic/', include('biostatistic.urls', namespace='biostatistic')),
                  path('archive/<int:year>/<int:month>/', ArchiveView.as_view(), name='archive'),
                  path('', generic.TemplateView.as_view(template_name="cover.html"), name='anasayfa'),
                  path('laboratuvarlar/', generic.TemplateView.as_view(template_name="laboratory.html"),
                       name='lab_home'),

                  path('laboratuvarlar/cbs/', generic.TemplateView.as_view(template_name="cover.html"), name="cbs"),

              ] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

urlpatterns += staticfiles_urlpatterns()
